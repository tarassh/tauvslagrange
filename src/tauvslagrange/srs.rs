use lambdaworks_math::{
    cyclic_group::IsGroup,
    elliptic_curve::{
        short_weierstrass::curves::bls12_381::{curve::BLS12381Curve, default_types::FrElement},
        traits::IsEllipticCurve,
    },
};
use rayon::prelude::*;

use crate::G1Point;

/// Generate SRS for a tau
pub fn generate_srs(n: usize, tau: FrElement) -> Vec<G1Point> {
    // Generate powers of tau: tau^1, tau^2, ..., tau^n
    let powers_of_tau = vandemonde_challenge(&tau, n - 1);

    let g1 = <BLS12381Curve as IsEllipticCurve>::generator();
    let mut tau_g1 = vec![g1; n];

    // Compute tau^i * g1 for i = 1, ..., n-1 in parallel
    tau_g1
        .par_iter_mut()
        .skip(1)
        .zip(&powers_of_tau)
        .for_each(|(g1, tau_i)| {
            *g1 = g1.operate_with_self(tau_i.representative());
        });

    tau_g1
}

/// Computes the powers of tau: tau^1, tau^2, ..., tau^n
fn vandemonde_challenge(x: &FrElement, n: usize) -> Vec<FrElement> {
    let mut powers = Vec::with_capacity(n);
    powers.push(x.clone());
    for i in 0..n - 1 {
        powers.push(x.pow(i as u64 + 2));
    }
    powers
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_vandemonde_challenge() {
        let challenge = vandemonde_challenge(&FrElement::from(2), 5);

        assert_eq!(
            challenge,
            vec![
                FrElement::from(2),
                FrElement::from(4),
                FrElement::from(8),
                FrElement::from(16),
                FrElement::from(32)
            ]
        );
    }
}
