# Memory in Hsp90 gating: discrete modes plus a slowly wandering closing rate

2026-10-15 · Roberto Covino (one-page summary of the main result after the independent audit; live version in Claude Docs)

Switching among a few discrete gating modes does not account for the memory in the open–closed dynamics of Hsp90; the closing rate must also wander, by about 20% over minutes.

## Data

Four apo yeast Hsp90 dimers were each followed for 6 h at 10 Hz with a plasmon ruler (Ye et al., Nano Lett. 2018; Zenodo 12697153). The signal has two levels, open and closed, separated by 11 standard deviations of the noise. The sequence of 25,571 open and closed dwells is therefore known without error, and since each frame is a snapshot, it contains all the information. Successive dwells are correlated (Spearman r = 0.35 between successive open times, decaying over 50–70 cycles), so the dynamics projected onto the open–closed coordinate are not Markovian (Vollmar et al., J. Phys. A 57, 365001, 2024).

## Open question

What is the hidden variable? It could be discrete: the molecule switches among a few gating modes, each with fixed opening and closing rates. Or the rates themselves could fluctuate continuously, as in Zwanzig's dynamic disorder. Both generate memory. The question is which of the two the data require.

## Approach

For any kinetic scheme with rate matrix Q, the probability of the observed frame sequence is a product of blocks of exp(QΔt), with Δt = 0.1 s. The likelihood is thus exact, including events shorter than a frame. The discrete reference is the best, by Bayesian evidence, of the ten discrete schemes we compared: five gating modes, ten states, 24 rate constants, detailed balance. We then let every closing rate depend on one coordinate z:

    k_OC^(m)(z) = k_OC^(m) exp(σz),   dz = −(z/τw) dt + (2/τw)^½ dW,   G_C(z) − G_O(z) = −σz

Here z diffuses in a harmonic well of unit variance with relaxation time τw, and the closed state shifts the well by σ (energies in kBT), which preserves detailed balance. On a grid, converged at 15 points, this is again a finite kinetic scheme, and the likelihood stays exact. We compared the models by maximum likelihood, by the evidence, and by predicting each molecule from the other three. The likelihoods and the evidence difference were reproduced independently, starting from the raw counts.

## Outcome

Adding these two parameters raises ln L by 77.5 and ln Z by 72. Molecules left out of the fit are predicted better by 76 in total. The fit gives σ = 0.22 (0.19–0.26) and τw = 300 s (200–450 s), 90% intervals. The same σ and τw, and gains of 77–94, appear on four different discrete skeletons, so the new term describes something the modes do not. A continuous z beats a two-state switch with as many or more parameters by 10 in ln L. The episodes that carry the gain show no change in intensity levels or noise.

The physics is simple. Since τw is about a hundred times the mean dwell, z is essentially frozen over some forty cycles. Two open times in the same mode are then correlated only because they share the same z:

    ρ = (e^{σ²} − 1)/(2e^{σ²} − 1) ≈ σ² ≈ 0.05

A correlation this weak is barely visible in dwell histograms or correlation functions. This is why it took the full likelihood of 25,000 dwells to detect it, and why our earlier simulation-based analysis on summary statistics had favored the discrete model.

Three findings limit the interpretation:

- **The effect is not uniform.** One molecule shows none (held-out Δln L = −7). Five of 48 half-hour segments carry 71% of the gain, more than in any of 40 data sets simulated from the fitted model, and the total gain is well below the model's own prediction (median 124).
- **The gating skeleton is too persistent.** The local closing rate loses its correlation within 6–8 min. Both models retain a correlation of 0.15 at 8 min, three standard deviations above the data.
- **Part of the signal is in the fastest events.** Removing all one-frame dwells lowers the gain from 78 to 48.

The conclusion is firm but limited. Fixed rates within a handful of discrete modes are ruled out, and the closing rate wanders by about 20% over minutes. But σ and τw are effective parameters, not yet properties of a physical coordinate. Four molecules cannot tell whether every molecule does this, or whether the wandering belongs to the protein or to its tether on gold. More molecules and the ATP data can.

Sources: numbers from `results/exact/` (fits, evidence, CV) and `results/audit/` (audit); details in NOTES.md, sections "Exact-likelihood model comparison" and "Independent audit of the main result".
