This is a from-scratch implementation of Binius, a polynomial commitment scheme built over binary tower fields instead of prime fields. Addition in these fields is XOR, so there are no carries or modular reduction to worry about, and multiplication recurses through the tower using a Karatsuba-style decomposition.

The code is meant to be read in this order:

binary_fields.py defines the field itself: element representation, addition, multiplication, and inversion.

utils.py builds the polynomial machinery on top of the field: multilinear evaluation, Lagrange interpolation, and the Reed-Solomon-style row extension.

merkle.py is a standalone Merkle tree, used to commit to and open the extended rows.

simple_binius.py ties the three together into the actual protocol: the prover commits to a polynomial and proves its evaluation at a point, and the verifier checks that proof.

driver.py is a guided walkthrough that exercises all of the above in order and prints what is happening at each step.

Run it with:

```
python driver.py
```

Everything here was built on top of the following material.

The two papers this follows directly:

- Diamond and Posen, "Succinct Arguments over Towers of Binary Fields" (the original Binius paper) — https://eprint.iacr.org/2023/1784
- Diamond and Posen, "Polylogarithmic Proofs for Multilinears over Binary Towers" (FRI-Binius) — https://eprint.iacr.org/2024/504

Background on sumcheck and multilinear extensions:

- Justin Thaler, *Proofs, Arguments, and Zero-Knowledge* — https://people.cs.georgetown.edu/jthaler/ProofsArgsAndZK.pdf
- Sum-Check Protocol and Multilinear Extensions, a shorter primer — https://risencrypto.github.io/Sumcheck/
- sumcheck.zksecurity.xyz, an interactive sumcheck/MLE tutorial in SageMath — https://sumcheck.zksecurity.xyz/pages/00-introduction.md
- ZK Whiteboard Sessions (ZK Hack), the ZKP and sumcheck modules

Talks:

- Ben Diamond, "Succinct Arguments over Towers of Binary Fields" — https://www.youtube.com/watch?v=eTCjVTWqjj0
- Jim Posen, "FRI-Binius: Polynomial Commitments for Tiny Binary Fields" — https://www.youtube.com/watch?v=0VcyaxA1gwc
- Kabir Peshawaria, "Introduction to Binius (State of the Art SNARK Proof System)" — https://www.youtube.com/watch?v=hB1STNfT920
- Kabir Peshawaria, "Introduction to Zerocheck (Part 1)" — https://www.youtube.com/watch?v=igz7MXsi2M0
- Kabir Peshawaria, "Introduction to Sumcheck (Part Two of Zerocheck Lecture)" — https://www.youtube.com/watch?v=hUJReqpjEvk
- Kabir Peshawaria, "Introduction to FRI Protocol" — https://www.youtube.com/watch?v=fGkvKXJbV_g
- Kabir Peshawaria, "Introduction to Additive NTT and Reed-Solomon Codes over Binary Fields" — https://www.youtube.com/watch?v=bAOJ2evWaek

Writeups:

- Vitalik Buterin, "Binius: highly efficient proofs over binary fields" — https://vitalik.eth.limo/general/2024/04/29/binius.html
- Vitalik Buterin, "STARKs, Part I: Proofs with Polynomials" and "Part II: Thank Goodness It's FRI-day" — https://vitalik.eth.limo/general/2017/11/09/starks_part_1.html, https://vitalik.eth.limo/general/2017/11/22/starks_part_2.html
- l2iterative, "Understanding Binius," Part I and Part II — https://hackmd.io/@l2iterative/binius, https://hackmd.io/@l2iterative/binius2
- Irreducible, "Binius: a Hardware-Optimized SNARK" — https://www.irreducible.com/posts/binius-hardware-optimized-snark
- Irreducible, "Binary Tower Fields are the Future of Verifiable Computing" — https://www.irreducible.com/posts/binary-tower-fields-are-the-future-of-verifiable-computing
- SECBIT, "Notes on FRI-Binius (Part I): Binary Towers" — https://secbit.io/blog/en/2024/10/31/binius-01/
- John D. Cook, "Why and how Bitcoin uses Merkle trees" — https://www.johndcook.com/blog/2025/10/28/bitcoin-merkle-trees/

Reference code, used as ground truth throughout:

- `ethereum/research/binius`, Vitalik's own simple/packed Binius implementation and the direct model for this build — https://github.com/ethereum/research/blob/master/binius/
- `IrreducibleOSS/binius-models`, the official minimal Python reference implementation — https://github.com/IrreducibleOSS/binius-models
- `IrreducibleOSS/binius`, the production Rust implementation, docs at docs.binius.xyz — https://github.com/IrreducibleOSS/binius

For more: awesome-binius, a curated index — https://github.com/kurtpan666/awesome-binius
