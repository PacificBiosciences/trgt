use super::Result;
use std::{
    collections::HashSet,
    env, fmt,
    path::{Path, PathBuf},
    str::FromStr,
    sync::{Arc, Mutex, OnceLock},
};
use url::Url;

// TODO: If the index is not local, it will automatically be downloaded to the local filesystem, we might want to have a post-processing step to delete the index after use...

fn redact_url(u: &Url) -> String {
    let mut redacted = u.clone();
    if !redacted.username().is_empty() || redacted.password().is_some() {
        let _ = redacted.set_username("");
        let _ = redacted.set_password(None);
    }
    if redacted.query().is_some() {
        redacted.set_query(None);
    }
    if redacted.fragment().is_some() {
        redacted.set_fragment(None);
    }
    redacted.to_string()
}

#[derive(Clone)]
pub enum Remote {
    Http(Url),
    S3(Url),
    Gcs(Url),
}

impl Remote {
    pub fn url(&self) -> &Url {
        match self {
            Remote::Http(u) | Remote::S3(u) | Remote::Gcs(u) => u,
        }
    }

    pub fn scheme_str(&self) -> &str {
        match self {
            Remote::Http(u) => u.scheme(), // "http" or "https"
            Remote::S3(_) => "s3",
            Remote::Gcs(_) => "gs",
        }
    }
}

impl fmt::Display for Remote {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", redact_url(self.url()))
    }
}

impl fmt::Debug for Remote {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Remote::Http(_) => f
                .debug_tuple("Http")
                .field(&redact_url(self.url()))
                .finish(),
            Remote::S3(_) => f.debug_tuple("S3").field(&redact_url(self.url())).finish(),
            Remote::Gcs(_) => f.debug_tuple("Gcs").field(&redact_url(self.url())).finish(),
        }
    }
}

#[derive(Clone)]
pub enum InputSource {
    Local(PathBuf),
    Remote(Remote),
}

#[derive(Clone, Debug)]
pub struct BamSources {
    sources: Vec<InputSource>,
    labels: Vec<Arc<str>>,
}

impl BamSources {
    pub fn new(primary: InputSource, supplemental: Option<InputSource>) -> Self {
        let sources = {
            let mut sources = Vec::with_capacity(1 + usize::from(supplemental.is_some()));
            sources.push(primary);
            if let Some(supplemental) = supplemental {
                sources.push(supplemental);
            }
            sources
        };
        let labels = sources
            .iter()
            .map(|source| Arc::<str>::from(source.to_string()))
            .collect();

        Self { sources, labels }
    }

    pub fn primary(&self) -> &InputSource {
        self.sources
            .first()
            .expect("BAM sources always contain a primary source")
    }

    pub fn as_slice(&self) -> &[InputSource] {
        &self.sources
    }

    pub fn iter_labeled(&self) -> impl Iterator<Item = (&InputSource, Arc<str>)> + '_ {
        self.sources.iter().zip(self.labels.iter().cloned())
    }
}

impl FromStr for InputSource {
    type Err = String;

    fn from_str(s: &str) -> Result<Self> {
        if let Ok(mut url) = Url::parse(s) {
            match url.scheme() {
                "file" => {
                    let p = url
                        .to_file_path()
                        .map_err(|_| "Bad file:// URL".to_string())?;
                    return if p.exists() {
                        Ok(InputSource::Local(p))
                    } else {
                        Err(format!("File does not exist: {}", p.display()))
                    };
                }
                "http" | "https" => {
                    return Ok(InputSource::Remote(Remote::Http(url)));
                }
                "s3" => {
                    return Ok(InputSource::Remote(Remote::S3(url)));
                }
                "gs" | "gcs" => {
                    if url.scheme() == "gcs" {
                        let _ = url.set_scheme("gs");
                    }
                    return Ok(InputSource::Remote(Remote::Gcs(url)));
                }
                _ => {}
            }
        }

        let p = Path::new(s);
        if p.exists() {
            Ok(InputSource::Local(p.to_path_buf()))
        } else {
            Err(format!("File does not exist: {}", p.display()))
        }
    }
}

impl fmt::Display for InputSource {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InputSource::Local(p) => write!(f, "{}", p.display()),
            InputSource::Remote(r) => write!(f, "{}", r),
        }
    }
}

impl fmt::Debug for InputSource {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            InputSource::Local(p) => f.debug_tuple("Local").field(p).finish(),
            InputSource::Remote(r) => f.debug_tuple("Remote").field(r).finish(),
        }
    }
}

static CA_ONCE: OnceLock<()> = OnceLock::new();
static PREFLIGHT_DONE: OnceLock<Mutex<HashSet<String>>> = OnceLock::new();

impl InputSource {
    pub fn file_stem(&self) -> Option<String> {
        match self {
            InputSource::Local(p) => p.file_stem()?.to_str().map(|s| s.to_owned()),
            InputSource::Remote(r) => {
                let u = r.url();
                let last = u.path_segments()?.last()?;
                if last.is_empty() {
                    return None;
                }
                let stem = last.rsplit_once('.').map(|(b, _)| b).unwrap_or(last);
                Some(stem.to_string())
            }
        }
    }

    /// Use this when wrapping errors from htslib that might leak URLs
    pub fn format_error(&self, context: &str, error: impl std::fmt::Display) -> String {
        let error_str = error.to_string();
        let redacted_error = if let InputSource::Remote(r) = self {
            let original_url = r.url().to_string();
            let redacted_url = redact_url(r.url());
            error_str.replace(&original_url, &redacted_url)
        } else {
            error_str
        };
        format!("{} {}: {}", context, self, redacted_error)
    }

    fn preflight_key(&self) -> Option<String> {
        match self {
            InputSource::Local(_) => None,
            InputSource::Remote(r) => {
                let scheme = r.scheme_str();
                let url = r.url();
                let host = url.host_str().unwrap_or_default().to_ascii_lowercase();
                let key = if let Some(port) = url.port() {
                    format!("{}://{}:{}", scheme, host, port)
                } else {
                    format!("{}://{}", scheme, host)
                };
                Some(key)
            }
        }
    }

    fn preflight_checks_inner(&self) -> Result<()> {
        if let InputSource::Remote(r) = self {
            match r {
                Remote::Http(u) if u.scheme() == "https" => {
                    CA_ONCE.get_or_init(set_system_ca_bundle);
                }
                Remote::S3(_) | Remote::Gcs(_) => {
                    CA_ONCE.get_or_init(set_system_ca_bundle);
                }
                _ => {}
            }

            let origin = self
                .preflight_key()
                .unwrap_or_else(|| r.scheme_str().to_string());

            match r {
                Remote::Http(u) => {
                    if u.scheme() == "https" && env::var_os("CURL_CA_BUNDLE").is_none() {
                        log::warn!(
                            "[TLS/CA] {origin}: No explicit CA bundle configured. \
                            HTSlib/libcurl will use the system's default trust store (native or file-based). \
                            If you hit TLS errors like 'SSL peer certificate ...', install system CA certificates \
                            or set CURL_CA_BUNDLE=/path/to/cert.pem."
                        );
                    }
                }
                Remote::S3(_) => {
                    let s3_vars = [
                        "AWS_ACCESS_KEY_ID",
                        "AWS_SECRET_ACCESS_KEY",
                        "AWS_SESSION_TOKEN",
                        "AWS_DEFAULT_REGION",
                        "AWS_PROFILE",
                        "AWS_SHARED_CREDENTIALS_FILE",
                    ];
                    if s3_vars.iter().all(|&v| env::var_os(v).is_none()) {
                        log::warn!("[S3 auth] {origin}: No AWS credential-related environment variables detected \
                             (AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY / AWS_SESSION_TOKEN / AWS_PROFILE / \
                              AWS_SHARED_CREDENTIALS_FILE / AWS_DEFAULT_REGION). \
                             Public buckets can usually be read anonymously; private buckets require credentials. \
                             Configure this through environment variables, shared config files (~/.aws/credentials and ~/.aws/config), or an instance/role/SSO provider. \
                             You may have to explicitly force anonymous reading, by setting AWS_SHARED_CREDENTIALS_FILE=/dev/null \
                             or use an equivalent https:// URL. TRGT does not change AWS_* automatically."
                        );
                    }
                }
                Remote::Gcs(_) => {
                    if env::var_os("GOOGLE_APPLICATION_CREDENTIALS").is_none() {
                        log::warn!("[GCS auth] {origin}: 'GOOGLE_APPLICATION_CREDENTIALS' is not set. \
                             Private buckets/objects require credentials, set GOOGLE_APPLICATION_CREDENTIALS=/path/to/key.json."
                        );
                    }
                }
            }
        }
        Ok(())
    }

    pub fn preflight_checks(&self) -> Result<()> {
        let Some(key) = self.preflight_key() else {
            return Ok(());
        };

        let set = PREFLIGHT_DONE.get_or_init(|| Mutex::new(HashSet::new()));
        {
            let guard = set.lock().unwrap_or_else(|p| p.into_inner());
            if guard.contains(&key) {
                return Ok(());
            }
        }

        let res = self.preflight_checks_inner();
        if res.is_ok() {
            let mut guard = set.lock().unwrap_or_else(|p| p.into_inner());
            guard.insert(key);
        }
        res
    }
}

fn set_system_ca_bundle() {
    if env::var_os("CURL_CA_BUNDLE").is_some()
        || env::var_os("TRGT_DISABLE_CA_AUTODETECT").is_some()
    {
        return;
    }

    const CANDIDATES: &[&str] = &[
        "/etc/ssl/certs/ca-certificates.crt", // Debian/Ubuntu/Gentoo
        "/etc/ssl/cert.pem",                  // Alpine/Arch/macOS/others
        "/opt/homebrew/etc/ca-certificates/cert.pem", // macOS Silicon
        "/usr/local/etc/ca-certificates/cert.pem", // macOS Intel
        "/usr/local/etc/openssl/cert.pem",    // macOS Legacy
        "/etc/ssl/ca-bundle.pem",             // OpenSUSE
        "/etc/pki/ca-trust/extracted/pem/tls-ca-bundle.pem", // CentOS/RHEL/Fedora
        "/etc/pki/tls/certs/ca-bundle.crt",   // Fedora/RHEL 6
        "/etc/pki/tls/cacert.pem",            // OpenELEC
    ];

    if let Some(found) = CANDIDATES.iter().find(|p| Path::new(p).exists()) {
        env::set_var("CURL_CA_BUNDLE", found);
        log::trace!("Using system CURL CA bundle found at {}", found);
        return;
    }

    log::debug!("No system CA bundle file found; relying on platform/libcurl defaults.");
}

pub trait PreflightExt {
    fn preflight_ext(&self) -> Result<()>;
}

impl PreflightExt for InputSource {
    #[inline]
    fn preflight_ext(&self) -> Result<()> {
        self.preflight_checks()
    }
}

impl PreflightExt for Option<InputSource> {
    #[inline]
    fn preflight_ext(&self) -> Result<()> {
        if let Some(inner) = self {
            inner.preflight_checks()?;
        }
        Ok(())
    }
}

#[macro_export]
macro_rules! preflight_fields {
    ($args:expr, $($field:ident),+ $(,)?) => {{
        use $crate::utils::input_source::PreflightExt as _;
        $( ($args).$field.preflight_ext()?; )+
        Ok::<(), std::string::String>(())
    }};
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::path::PathBuf;
    use url::Url;

    fn with_basic_auth(mut url: Url, user: &str, pass: &str) -> Url {
        url.set_username(user).unwrap();
        url.set_password(Some(pass)).unwrap();
        url
    }
    fn with_username_only(mut url: Url, user: &str) -> Url {
        url.set_username(user).unwrap();
        url
    }
    fn with_query(mut url: Url, k: &str, v: &str) -> Url {
        url.query_pairs_mut().append_pair(k, v);
        url
    }
    fn with_fragment(mut url: Url, frag: &str) -> Url {
        url.set_fragment(Some(frag));
        url
    }

    #[test]
    fn test_redact_url_basic_auth() {
        let base = Url::parse("s3://bucket/path/sample.bam").unwrap();
        let url = with_basic_auth(base, "myuser", "secretkey");
        let redacted = redact_url(&url);
        assert_eq!(redacted, "s3://bucket/path/sample.bam");
        assert!(!redacted.contains("myuser"));
        assert!(!redacted.contains("secretkey"));
        assert!(!redacted.contains('@'));
    }

    #[test]
    fn test_redact_url_query_string() {
        let base = Url::parse("https://example.com/sample.bam").unwrap();
        let url = with_query(base, "token", "secret123");
        let redacted = redact_url(&url);
        assert_eq!(redacted, "https://example.com/sample.bam");
        assert!(!redacted.contains('?'));
        assert!(!redacted.contains("token"));
        assert!(!redacted.contains("secret123"));
    }

    #[test]
    fn test_redact_url_fragment() {
        let base = Url::parse("https://example.com/sample.bam").unwrap();
        let url = with_fragment(base, "sensitive-fragment");
        let redacted = redact_url(&url);
        assert_eq!(redacted, "https://example.com/sample.bam");
        assert!(!redacted.contains('#'));
        assert!(!redacted.contains("sensitive-fragment"));
    }

    #[test]
    fn test_redact_url_no_credentials() {
        let url = Url::parse("s3://bucket/path/sample.bam").unwrap();
        let redacted = redact_url(&url);
        assert_eq!(redacted, "s3://bucket/path/sample.bam");
    }

    #[test]
    fn test_redact_url_only_username() {
        let base = Url::parse("https://example.com/sample.bam").unwrap();
        let url = with_username_only(base, "user");
        let redacted = redact_url(&url);
        assert_eq!(redacted, "https://example.com/sample.bam");
        assert!(!redacted.contains("user"));
        assert!(!redacted.contains('@'));
    }

    #[test]
    fn bam_sources_preserve_primary_and_supplemental_order() {
        let primary = InputSource::Local(PathBuf::from("/inputs/hifi.bam"));
        let supplemental = InputSource::Local(PathBuf::from("/inputs/fail.bam"));

        let sources = BamSources::new(primary.clone(), Some(supplemental.clone()));

        assert_eq!(sources.as_slice().len(), 2);
        assert_eq!(sources.primary().to_string(), primary.to_string());
        assert_eq!(sources.as_slice()[1].to_string(), supplemental.to_string());

        let labels = sources
            .iter_labeled()
            .map(|(_, label)| label.to_string())
            .collect::<Vec<_>>();
        assert_eq!(labels, vec![primary.to_string(), supplemental.to_string()]);
    }

    #[test]
    fn test_format_error_redacts_url() {
        let base = Url::parse("s3://bucket/sample.bam").unwrap();
        let url = with_basic_auth(base, "myuser", "secretpass");
        let src = InputSource::Remote(Remote::S3(url.clone()));
        let error = format!("Failed to open {}: Permission denied", url);
        let result = src.format_error("Failed to open", &error);
        assert!(result.contains("s3://bucket/sample.bam"));
        assert!(!result.contains("myuser"));
        assert!(!result.contains("secretpass"));
        assert!(!result.contains('@'));
    }

    #[test]
    fn test_format_error_local_unchanged() {
        let src = InputSource::Local(PathBuf::from("/path/to/file.bam"));
        let remote = with_basic_auth(Url::parse("s3://bucket/file").unwrap(), "user", "pass");
        let error = format!("Some error message with {}", remote);
        let result = src.format_error("Failed", &error);
        let mut prefix = remote.clone();
        prefix.set_path("");
        prefix.set_query(None);
        prefix.set_fragment(None);
        let expect_prefix = prefix.as_str().trim_end_matches('/').to_string();
        assert!(result.contains(&expect_prefix));
    }

    #[test]
    fn test_format_error_with_query() {
        let base = Url::parse("https://example.com/sample.bam").unwrap();
        let url = with_query(base, "apikey", "secret123");
        let src = InputSource::Remote(Remote::Http(url.clone()));
        let error = format!("Failed to fetch {}: timeout", url);
        let result = src.format_error("Failed to fetch", &error);
        assert!(result.contains("https://example.com/sample.bam"));
        assert!(!result.contains("apikey"));
        assert!(!result.contains("secret123"));
    }

    #[test]
    fn test_redact_url_combined() {
        let base = Url::parse("https://example.com/sample.bam").unwrap();
        let url = with_fragment(
            with_query(with_basic_auth(base, "user", "pass"), "token", "abc"),
            "section",
        );
        let redacted = redact_url(&url);
        assert_eq!(redacted, "https://example.com/sample.bam");
        assert!(!redacted.contains("user"));
        assert!(!redacted.contains("pass"));
        assert!(!redacted.contains('@'));
        assert!(!redacted.contains('?'));
        assert!(!redacted.contains('#'));
        assert!(!redacted.contains("token"));
        assert!(!redacted.contains("abc"));
        assert!(!redacted.contains("section"));
    }
}
