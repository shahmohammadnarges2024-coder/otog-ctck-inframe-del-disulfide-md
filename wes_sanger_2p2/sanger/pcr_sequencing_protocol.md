# PCR and Sequencing Protocol for OTOG Exon 55

**PCR reaction (25 µL):**
- 1× PCR buffer
- 1.5 mM MgCl₂
- 0.2 mM each dNTP
- 0.2 µM each primer
- 1 U Taq DNA polymerase
- 20–50 ng genomic DNA

**Touchdown cycling program (to accommodate ~10 °C primer Tm difference):**
- 94 °C for 4 min
- 14 touchdown cycles: 94 °C for 40 s; annealing 64 °C → 57 °C (−0.5 °C/cycle) for 30 s; 72 °C for 40 s
- 18 standard cycles: 94 °C for 40 s; 57 °C for 30 s; 72 °C for 40 s
- Final extension: 72 °C for 10 min

**Amplicon check:**
- 2% agarose gel with SYBR™ Safe
- Size: ~385 bp (GeneRuler 100 bp DNA Ladder, Thermo Scientific)

**Sequencing:**
- Purification: ExoSAP-IT™ (Applied Biosystems)
- Cycle sequencing: BigDye™ Terminator v3.1, PCR primers
- Capillary electrophoresis: ABI 3500 Genetic Analyzer
- Analysis: CodonCode Aligner v11.0.2, aligned to GRCh38
