use lambdaworks_math::{
    elliptic_curve::{
        short_weierstrass::curves::bls12_381::curve::BLS12381Curve, traits::IsEllipticCurve,
    },
    unsigned_integer::element::UnsignedInteger,
};
use serde::{Deserialize, Serialize};

use crate::G1Point;

#[derive(Debug, Serialize, Deserialize)]
pub struct SerializedSRS {
    pub points: Vec<(String, String)>,
}

impl From<Vec<G1Point>> for SerializedSRS {
    fn from(srs: Vec<G1Point>) -> Self {
        let affined = srs
            .iter()
            .map(|p| (p.to_affine().x().to_string(), p.to_affine().y().to_string()))
            .collect::<Vec<_>>();

        SerializedSRS { points: affined }
    }
}

impl SerializedSRS {
    pub fn to_ec_points(self) -> Vec<G1Point> {
        self.points
            .iter()
            .map(|(x, y)| {
                let x = UnsignedInteger::from_hex_unchecked(x);
                let y = UnsignedInteger::from_hex_unchecked(y);

                <BLS12381Curve as IsEllipticCurve>::create_point_from_affine(
                    (&x).into(),
                    (&y).into(),
                )
                .unwrap()
            })
            .collect()
    }
}

impl SerializedSRS {
    pub fn dump(&self, file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
        let serialized_data = serde_json::to_string(&self.points)?;
        std::fs::write(file_path, serialized_data)?;

        Ok(())
    }

    pub fn load(file_path: &str) -> Result<Self, Box<dyn std::error::Error>> {
        let serialized_data = std::fs::read_to_string(file_path)?;
        let points: Vec<(String, String)> = serde_json::from_str(&serialized_data)?;

        Ok(SerializedSRS { points })
    }
}
