# Next external expression unit: admission and frozen execution boundary

This contract does not admit a new source or claim validation. It follows the existing negative external-source audit and the already-admitted GSE63789 count-ratio transport.

## Construct intervention source (Nature16509 candidate)
Require the actual construct sequence, unique construct ID, assay units and scoring direction, host strain, promoter/vector and replicate identifiers. Separate natural genes from synonymous variants. Require a mapping from construct to template protein to prevent sequence-neighbor leakage. The currently available predictors cannot automatically interpret promoter/T7 context or an ordinal expression score as endogenous protein abundance. If the source has multiple variants per template, hold out whole template groups rather than randomly splitting constructs. Source admission must establish license/provenance and sequences/outcomes, not just paper claims. Freeze all endpoint transformations and source-only model choices before outcome scoring. Paired synonymous-effect validation requires parent-variant links and consistent assay context; natural-gene cross-source correlation is not a substitute.

## Conditional Ribo/RNA source (GSE182100 candidate)
Require actual accession/run-to-condition/replicate mappings, feature identifiers, count or normalized endpoint definitions and host strain/reference assembly. Keep repeated genes across conditions marked as repeated identities. Confirm RNA and footprint libraries are paired and do not mix a condition-specific ratio with absolute protein yield. Prespecify eligible conditions and all exclusions before target scoring, retaining all conditions rather than selecting a favorable one. Same-gene condition transport is separate from held-out-gene generalization. Biological-replicate uncertainty needs actual paired replicates; a single ratio and gene resampling cannot certify it.

## Work possible before admission
Confirm feature geometry and predictor input requirements from current code; prepare parsers only against documented schemas, without inventing a file format. No target fitting, calibration, construct optimization or design-benefit assertion is permitted by this contract. If no usable sequence/endpoint pairing exists, return the precise gap instead of translating the source into the nearest old endpoint.

## Completion boundary
Report source URLs, file hashes, row/identity units, overlap and selection history, exact frozen model/metric scope, all negatives and remaining gates. No new source is assumed independent merely because it has a new accession or paper.
