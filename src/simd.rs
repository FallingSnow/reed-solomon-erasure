// https://rust-lang.github.io/rfcs/2045-target-feature.html

use core::arch::x86_64::{__m128i, __m256i};

static mut VECTOR_SIZE: usize = 16;

union V {
    m128i: __m128i,
    m256i: __m256i,
}

impl From<__m128i> for V {
    fn from(value: __m128i) -> Self {
        Self { m128i: value }
    }
}

impl From<__m256i> for V {
    fn from(value: __m256i) -> Self {
        Self { m256i: value }
    }
}

/** And */
#[inline]
unsafe fn and_rust(a: &V, b: &V) -> V {
    todo!()
}

#[target_feature(enable = "sse2")]
#[inline]
unsafe fn and_sse2(a: &V, b: &V) -> V {
    core::arch::x86_64::_mm_and_si128(a.m128i, b.m128i).into()
}

#[target_feature(enable = "avx2")]
#[inline]
unsafe fn and_avx2(a: &V, b: &V) -> V {
    core::arch::x86_64::_mm256_and_si256(a.m256i, b.m256i).into()
}

static mut and_v: unsafe fn(a: &V, b: &V) -> V = and_rust;

/** Xor */
#[inline]
unsafe fn xor_rust(a: &V, b: &V) -> V {
    todo!()
}

#[target_feature(enable = "sse2")]
#[inline]
unsafe fn xor_sse2(a: &V, b: &V) -> V {
    core::arch::x86_64::_mm_xor_si128(a.m128i, b.m128i).into()
}

#[target_feature(enable = "avx2")]
#[inline]
unsafe fn xor_avx2(a: &V, b: &V) -> V {
    core::arch::x86_64::_mm256_xor_si256(a.m256i, b.m256i).into()
}

static mut xor_v: unsafe fn(a: &V, b: &V) -> V = xor_rust;

/** Store */
#[inline]
unsafe fn store_rust(out: &mut [u8], vec: &V) {
    todo!()
}

#[target_feature(enable = "sse2")]
#[inline]
unsafe fn store_sse2(out: &mut [u8], vec: &V) {
    core::arch::x86_64::_mm_storeu_si128(out.as_mut_ptr() as *mut __m128i, vec.m128i).into()
}

#[target_feature(enable = "avx2")]
#[inline]
unsafe fn store_avx2(out: &mut [u8], vec: &V) {
    core::arch::x86_64::_mm256_storeu_si256(out.as_mut_ptr() as *mut __m256i, vec.m256i).into()
}

static mut store_v: unsafe fn(out: &mut [u8], vec: &V) = store_rust;

/** Load */

#[target_feature(enable = "avx2")]
#[inline]
unsafe fn load_128_avx2(input: &[u8]) -> V {
    let r = load_sse2(input);
    core::arch::x86_64::_mm256_broadcastsi128_si256(r.m128i).into()
}

static mut load_128_v: unsafe fn(input: &[u8]) -> V = load_rust;

/** Load */
#[inline]
unsafe fn load_rust(input: &[u8]) -> V {
    todo!()
}

#[target_feature(enable = "sse2")]
#[inline]
unsafe fn load_sse2(input: &[u8]) -> V {
    core::arch::x86_64::_mm_loadu_si128(input.as_ptr() as *const __m128i).into()
}

#[target_feature(enable = "avx2")]
#[inline]
unsafe fn load_avx2(input: &[u8]) -> V {
    core::arch::x86_64::_mm256_loadu_si256(input.as_ptr() as *const __m256i).into()
}

static mut load_v: unsafe fn(input: &[u8]) -> V = load_rust;

/** Splat */
#[inline]
unsafe fn splat_rust(input: u8) -> V {
    todo!()
}

#[target_feature(enable = "sse2")]
#[inline]
unsafe fn splat_sse2(input: u8) -> V {
    core::arch::x86_64::_mm_set1_epi8(input as i8).into()
}

#[target_feature(enable = "avx2")]
#[inline]
unsafe fn splat_avx2(input: u8) -> V {
    core::arch::x86_64::_mm256_set1_epi8(input as i8).into()
}

static mut splat_v: unsafe fn(input: u8) -> V = splat_rust;

/** Shift */
#[inline]
unsafe fn shift_rust(vec: &V) -> V {
    todo!()
}

#[target_feature(enable = "sse2")]
#[inline]
unsafe fn shift_sse2(vec: &V) -> V {
    core::arch::x86_64::_mm_srli_epi64(vec.m128i, 4).into()
}

#[target_feature(enable = "avx2")]
#[inline]
unsafe fn shift_avx2(vec: &V) -> V {
    core::arch::x86_64::_mm256_srli_epi64(vec.m256i, 4).into()
}

static mut shift_v: unsafe fn(vec: &V) -> V = shift_rust;

/** Shuffle */
#[inline]
unsafe fn shuffle_rust(a: &V, b: &V) -> V {
    todo!()
}

#[target_feature(enable = "ssse3")]
#[inline]
unsafe fn shuffle_ssse3(a: &V, b: &V) -> V {
    core::arch::x86_64::_mm_shuffle_epi8(a.m128i, b.m128i).into()
}

#[target_feature(enable = "avx2")]
#[inline]
unsafe fn shuffle_avx2(a: &V, b: &V) -> V {
    core::arch::x86_64::_mm256_shuffle_epi8(a.m256i, b.m256i).into()
}

static mut shuffle_v: unsafe fn(a: &V, b: &V) -> V = shuffle_rust;

#[derive(PartialEq, Eq, PartialOrd, Ord)]
pub enum InstructionSet {
    SSE2,
    SSSE3,
    AVX2,
    AVX512,
}

unsafe fn set_sse2() {
    VECTOR_SIZE = 16;
    and_v = and_sse2;
    xor_v = xor_sse2;
    store_v = store_sse2;
    load_128_v = load_sse2;
    load_v = load_sse2;
    splat_v = splat_sse2;
    shift_v = shift_sse2;
}
unsafe fn set_ssse3() {
    shuffle_v = shuffle_ssse3;
}
unsafe fn set_avx2() {
    VECTOR_SIZE = 32;
    and_v = and_avx2;
    xor_v = xor_avx2;
    store_v = store_avx2;
    load_128_v = load_128_avx2;
    load_v = load_avx2;
    splat_v = splat_avx2;
    shift_v = shift_avx2;
    shuffle_v = shuffle_avx2;
}

pub fn set_instruction_set(simd: Option<InstructionSet>) {
    if let Some(set) = simd {
        use InstructionSet::*;
        unsafe {
            match set {
                SSE2 => set_sse2(),
                SSSE3 => {
                    set_sse2();
                    set_ssse3();
                }
                AVX2 => set_avx2(),
                AVX512 => unimplemented!(),
            }
        }
        return;
    }

    // run-time feature detection:
    #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
    unsafe {
        if is_x86_feature_detected!("avx2") {
            // println!("Using AVX2");
            set_avx2();
        } else if is_x86_feature_detected!("sse2") {
            // println!("Using SSE2");
            set_sse2();
            if is_x86_feature_detected!("ssse3") {
                set_ssse3();
            }
        }
    }

    // NEON
    #[cfg(any(target_arch = "arm", target_arch = "aarch64"))]
    {}
}

pub unsafe fn reedsolomon_gal_mul(
    low: &[u8; 16],
    high: &[u8; 16],
    input: &[u8],
    output: &mut [u8],
    xor: bool,
) -> usize {
    let low_mask_unpacked = splat_v(0x0f);
    let low_vector = load_128_v(low);
    let high_vector = load_128_v(high);

    let mut done: usize = 0;
    for _ in 0..(input.len() / VECTOR_SIZE) {
        let in_x = load_v(&input[done..]);
        let old = load_v(&output[done..]);

        let result = {
            let low_input = and_v(&in_x, &low_mask_unpacked);
            let in_x_shifted = shift_v(&in_x);
            let high_input = and_v(&in_x_shifted, &low_mask_unpacked);

            let mul_low_part = shuffle_v(&low_vector, &low_input);
            let mul_high_part = shuffle_v(&high_vector, &high_input);

            let mut new = xor_v(&mul_low_part, &mul_high_part);
            if xor {
                new = xor_v(&new, &old);
            }
            new
        };

        store_v(&mut output[done..], &result);

        done += VECTOR_SIZE;
    }

    // dbg!(done);

    done
}
