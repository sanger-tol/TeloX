fn is_dinucleotide_repeat_bytes(seq: &[u8]) -> bool {
    let len = seq.len();
    if len < 4 || len % 2 != 0 {
        return false;
    }
    let first = seq[0];
    let second = seq[1];
    if first == second {
        return false; // not a dinucleotide repeat if both are the same
    }
    for i in (0..len).step_by(2) {
        if seq[i] != first || seq[i + 1] != second {
            return false;
        }
    }
    true
}

fn main() {
    let test_cases = [
        ("ACACAC", true),
        ("ATATATAT", true),
        ("CGCGCG", true),
        ("AACCAACCAA", false), // 4-mer repeat, not dinucleotide
        ("ABCDEF", false),
        ("AAAA", false), // homopolymer
        ("ACACACG", false), // not perfect repeat
    ];
    
    for (seq, expected) in &test_cases {
        let result = is_dinucleotide_repeat_bytes(seq.as_bytes());
        println!("{}: {} (expected {})", seq, result, expected);
        if result != *expected {
            println!("  ❌ MISMATCH!");
        } else {
            println!("  ✅ OK");
        }
    }
}
