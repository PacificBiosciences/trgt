use std::fmt::{Binary, Display};

pub type Result<T> = std::result::Result<T, String>;

/// Returns the FNV-1a hash of a run of bytes.
///
/// It is written out rather than taken from the standard library because the values it gives decide
/// what a run writes and what a later run compares against, so they have to stay the same whatever
/// version of Rust builds TRGT. The values are therefore part of what TRGT's output files mean, and
/// the test below pins them.
pub fn fnv1a64(bytes: &[u8]) -> u64 {
    let mut hash: u64 = 0xcbf2_9ce4_8422_2325;
    for byte in bytes {
        hash ^= u64::from(*byte);
        hash = hash.wrapping_mul(0x0000_0100_0000_01b3);
    }
    hash
}

pub fn handle_error_and_exit(err: String) -> ! {
    log::error!("{}", err);
    std::process::exit(1);
}

pub fn format_number_with_commas<T>(n: T) -> String
where
    T: Display + Binary,
{
    let s = n.to_string();
    let (sign, digits) = s.strip_prefix('-').map_or(("", s.as_str()), |d| ("-", d));

    if let 0..=3 = digits.len() {
        return s;
    }

    let mut result = String::with_capacity(digits.len() + (digits.len() - 1) / 3 + sign.len());
    for (digit_count, c) in digits.chars().rev().enumerate() {
        if digit_count > 0 && digit_count % 3 == 0 {
            result.push(',');
        }
        result.push(c);
    }

    result = result.chars().rev().collect();
    if !sign.is_empty() {
        result.insert_str(0, sign);
    }

    result
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_number_with_commas_i8() {
        assert_eq!(format_number_with_commas(0i8), "0");
        assert_eq!(format_number_with_commas(100i8), "100");
        assert_eq!(format_number_with_commas(-100i8), "-100");
        assert_eq!(format_number_with_commas(i8::MAX), "127");
        assert_eq!(format_number_with_commas(i8::MIN), "-128");
    }

    #[test]
    fn test_format_number_with_commas_u8() {
        assert_eq!(format_number_with_commas(0u8), "0");
        assert_eq!(format_number_with_commas(100u8), "100");
        assert_eq!(format_number_with_commas(u8::MAX), "255");
    }

    #[test]
    fn test_format_number_with_commas_i16() {
        assert_eq!(format_number_with_commas(0i16), "0");
        assert_eq!(format_number_with_commas(1_000i16), "1,000");
        assert_eq!(format_number_with_commas(-1_000i16), "-1,000");
        assert_eq!(format_number_with_commas(i16::MAX), "32,767");
        assert_eq!(format_number_with_commas(i16::MIN), "-32,768");
    }

    #[test]
    fn test_format_number_with_commas_u16() {
        assert_eq!(format_number_with_commas(0u16), "0");
        assert_eq!(format_number_with_commas(1_000u16), "1,000");
        assert_eq!(format_number_with_commas(u16::MAX), "65,535");
    }

    #[test]
    fn test_format_number_with_commas_i32() {
        assert_eq!(format_number_with_commas(0i32), "0");
        assert_eq!(format_number_with_commas(10_000i32), "10,000");
        assert_eq!(format_number_with_commas(-10_000i32), "-10,000");
        assert_eq!(format_number_with_commas(i32::MAX), "2,147,483,647");
        assert_eq!(format_number_with_commas(i32::MIN), "-2,147,483,648");
    }

    #[test]
    fn test_format_number_with_commas_u32() {
        assert_eq!(format_number_with_commas(0u32), "0");
        assert_eq!(format_number_with_commas(10_000u32), "10,000");
        assert_eq!(format_number_with_commas(u32::MAX), "4,294,967,295");
    }

    #[test]
    fn test_format_number_with_commas_i64() {
        assert_eq!(format_number_with_commas(0i64), "0");
        assert_eq!(format_number_with_commas(1_000_000i64), "1,000,000");
        assert_eq!(format_number_with_commas(-1_000_000i64), "-1,000,000");
        assert_eq!(
            format_number_with_commas(i64::MAX),
            "9,223,372,036,854,775,807"
        );
        assert_eq!(
            format_number_with_commas(i64::MIN),
            "-9,223,372,036,854,775,808"
        );
    }

    #[test]
    fn test_format_number_with_commas_u64() {
        assert_eq!(format_number_with_commas(0u64), "0");
        assert_eq!(format_number_with_commas(1_000_000u64), "1,000,000");
        assert_eq!(
            format_number_with_commas(u64::MAX),
            "18,446,744,073,709,551,615"
        );
    }

    #[test]
    fn test_format_number_with_commas_usize() {
        assert_eq!(format_number_with_commas(0usize), "0");
        assert_eq!(
            format_number_with_commas(1_234_567_890usize),
            "1,234,567,890"
        );
        assert_eq!(
            format_number_with_commas(usize::MAX),
            "18,446,744,073,709,551,615"
        );
    }

    #[test]
    fn test_format_number_with_commas_isize() {
        assert_eq!(format_number_with_commas(0isize), "0");
        assert_eq!(
            format_number_with_commas(-1_234_567_890isize),
            "-1,234,567,890"
        );
        assert_eq!(
            format_number_with_commas(isize::MAX),
            "9,223,372,036,854,775,807"
        );
        assert_eq!(
            format_number_with_commas(isize::MIN),
            "-9,223,372,036,854,775,808"
        );
    }

    #[test]
    fn fnv1a64_gives_the_values_trgt_output_files_are_written_with() {
        // Changing these values would make TRGT refuse to resume every output file written before
        // the change, and would hand every deep locus a different sample of its reads.
        assert_eq!(fnv1a64(b""), 0xcbf2_9ce4_8422_2325);
        assert_eq!(fnv1a64(b"a"), 0xaf63_dc4c_8601_ec8c);
        assert_eq!(fnv1a64(b"foobar"), 0x8594_4171_f739_67e8);
    }

    #[test]
    fn fnv1a64_tells_different_bytes_apart() {
        assert_ne!(fnv1a64(b"chr1:100-200:TR1"), fnv1a64(b"chr1:100-200:TR2"));
        assert_ne!(fnv1a64(b"chr1:100-200:TR1"), fnv1a64(b"chr1:100-201:TR1"));
        assert_eq!(fnv1a64(b"chr1:100-200:TR1"), fnv1a64(b"chr1:100-200:TR1"));
    }
}
