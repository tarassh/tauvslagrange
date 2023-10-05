# Powers of Tau vs Lagrange Basis

## Results

```
*******************************
*                             *
*  Powers of Tau vs Lagrange  *
*                             *
*******************************


Options:
1. Run commitment with pre-generated SRS
2. Generate new SRS
3. Exit
> 1


------------ Setup ------------
Loading powers of tau ...
Loading powers of tau - Elapsed: 99.25125ms
Loading powers of tau in Lagrange basis ...
Loading powers of tau in Lagrange basis - Elapsed: 65.753583ms
Polynomial Generation ...
Polynomial Generation - Elapsed: 2.853042ms


------------ Prover ------------
Witness Generation ...
Witness Generation - Elapsed: 2.573541ms
Commitment Calculation (Powers of Tau) ...
Commitment Calculation (Powers of Tau) - Elapsed: 1.000773667s
Commitment Calculation (Lagrange) ...
Commitment Calculation (Lagrange) - Elapsed: 943.58275ms


------------ Result ------------
Commitment[t] G1: (0x967dbd96321619d880a2faa30118ce2466362e14d850ea45151479da9deecb9dfe42451ff1f7a6324fb3ad6a89fe026, 0x14f969cd42ed59a9b911849e90bbe8542c6555f9682860d2920b1a5e29e5ee1387ead3e7c0b58778219539d31ca33b67)
Commitment[l] G1: (0x967dbd96321619d880a2faa30118ce2466362e14d850ea45151479da9deecb9dfe42451ff1f7a6324fb3ad6a89fe026, 0x14f969cd42ed59a9b911849e90bbe8542c6555f9682860d2920b1a5e29e5ee1387ead3e7c0b58778219539d31ca33b67)
```

## Build

```
cargo build --release
```

## Run

```
cargo run --release
```

