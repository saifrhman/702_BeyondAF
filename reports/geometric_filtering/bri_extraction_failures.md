# PDB707K BRI extraction failures

## Run summary

- Production manifest: 461,165 model-1 chains
- Successful BRI vectors: 460,920
- Failed chains: 245
- Success rate: 99.9469%
- SLURM array tasks completed: 210
- Extraction command: `pdb_assess pdb707k-extract`

## Failure categories

| Failure category | Chains | Unique PDB entries |
| --- | ---: | ---: |
| Residue-count mismatch | 218 | 46 |
| Invariant not generated | 27 | 15 |
| Other | 0 | 0 |

## Interpretation

The residue-count mismatch cases occur when the current mmCIF coordinates
produce fewer invariant residues than expected by the historical PDB707K
manifest.

The invariant-not-generated cases require individual inspection. Some appear
to involve differences between historical chain identifiers and the current
label/auth chain mapping.

No coordinates were padded, fabricated or silently shortened. These chains
will be excluded unless a correction can be validated unambiguously.

## Exact failed-chain list

| pdb_id | model_id | chain_id | start_residue | chain_length | failure_type | error |
| --- | --- | --- | --- | --- | --- | --- |
| 1AAM | 1 | A | 1 | 396 | residue_count_mismatch | ValueError('Expected 396 invariant residues, got 395') |
| 1CRX | 1 | F | 1 | 322 | residue_count_mismatch | ValueError('Expected 322 invariant residues, got 321') |
| 1J0E | 1 | A | 1 | 341 | residue_count_mismatch | ValueError('Expected 341 invariant residues, got 340') |
| 1J0E | 1 | B | 1 | 341 | residue_count_mismatch | ValueError('Expected 341 invariant residues, got 340') |
| 1J0E | 1 | C | 1 | 341 | residue_count_mismatch | ValueError('Expected 341 invariant residues, got 340') |
| 1J0E | 1 | D | 1 | 341 | residue_count_mismatch | ValueError('Expected 341 invariant residues, got 340') |
| 1SFT | 1 | A | 2 | 382 | residue_count_mismatch | ValueError('Expected 382 invariant residues, got 381') |
| 1SFT | 1 | B | 2 | 380 | residue_count_mismatch | ValueError('Expected 380 invariant residues, got 379') |
| 2AAA | 1 | A | 1 | 476 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 2NBI | 1 | A | 1 | 494 | residue_count_mismatch | ValueError('Expected 494 invariant residues, got 485') |
| 2OMG | 1 | A | 1 | 21 | residue_count_mismatch | ValueError('Expected 21 invariant residues, got 20') |
| 2OMG | 1 | D | 1 | 28 | residue_count_mismatch | ValueError('Expected 28 invariant residues, got 27') |
| 2OMG | 1 | E | 1 | 21 | residue_count_mismatch | ValueError('Expected 21 invariant residues, got 20') |
| 2OMG | 1 | F | 1 | 28 | residue_count_mismatch | ValueError('Expected 28 invariant residues, got 27') |
| 2Y1L | 1 | H | 1 | 4 | residue_count_mismatch | ValueError('Expected 4 invariant residues, got 3') |
| 2Z1Y | 1 | B | 4 | 394 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 4TKY | 1 | E | 1 | 7 | residue_count_mismatch | ValueError('Expected 7 invariant residues, got 6') |
| 4TKY | 1 | F | 1 | 7 | residue_count_mismatch | ValueError('Expected 7 invariant residues, got 6') |
| 4TKY | 1 | G | 1 | 7 | residue_count_mismatch | ValueError('Expected 7 invariant residues, got 6') |
| 4TKY | 1 | H | 1 | 7 | residue_count_mismatch | ValueError('Expected 7 invariant residues, got 6') |
| 5CGH | 1 | CA | 4 | 1 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 5CGH | 1 | DA | 4 | 1 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 5LF6 | 1 | K | 1 | 2 | residue_count_mismatch | ValueError('Expected 2 invariant residues, got 1') |
| 5LF6 | 1 | Z | 1 | 2 | residue_count_mismatch | ValueError('Expected 2 invariant residues, got 1') |
| 5TGM | 1 | M | 2 | 155 | residue_count_mismatch | ValueError('Expected 155 invariant residues, got 154') |
| 5TGM | 1 | N | 20 | 124 | residue_count_mismatch | ValueError('Expected 124 invariant residues, got 105') |
| 5TGM | 1 | NC | 2 | 146 | residue_count_mismatch | ValueError('Expected 146 invariant residues, got 145') |
| 5TGM | 1 | OC | 20 | 124 | residue_count_mismatch | ValueError('Expected 124 invariant residues, got 105') |
| 5TGM | 1 | Q | 8 | 124 | residue_count_mismatch | ValueError('Expected 124 invariant residues, got 121') |
| 5TGM | 1 | RC | 4 | 135 | residue_count_mismatch | ValueError('Expected 135 invariant residues, got 132') |
| 5TGM | 1 | TD | 24 | 231 | residue_count_mismatch | ValueError('Expected 231 invariant residues, got 208') |
| 5U8X | 1 | A | 30 | 497 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 5U8X | 1 | C | 27 | 500 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 5U8X | 1 | D | 30 | 497 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 6D6V | 1 | B | 511 | 185 | residue_count_mismatch | ValueError('Expected 185 invariant residues, got 145') |
| 6NSR | 1 | E | 114 | 1 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 6QNX | 1 | B | 321 | 74 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 6QNX | 1 | C | 223 | 9 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 6TVW | 1 | B | 7 | 9 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 6XND | 1 | A | 3 | 121 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 6XND | 1 | B | 3 | 121 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 6XND | 1 | C | 3 | 121 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 6XND | 1 | D | 3 | 121 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 7CLR | 1 |  | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | AA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | BA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | CA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | DA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | EA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | FA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | GA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | HA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | IA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | JA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | KA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | LA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | MA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | OA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | PA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | QA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | RA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | SA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | TA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | UA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | VA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | WA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | XA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | YA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7CLR | 1 | ZA | 22 | 211 | residue_count_mismatch | ValueError('Expected 211 invariant residues, got 200') |
| 7ED1 | 1 | A | 2 | 223 | residue_count_mismatch | ValueError('Expected 223 invariant residues, got 222') |
| 7ED1 | 1 | B | 2 | 223 | residue_count_mismatch | ValueError('Expected 223 invariant residues, got 222') |
| 7K0I | 1 | B | 53 | 492 | residue_count_mismatch | ValueError('Expected 492 invariant residues, got 491') |
| 7K0I | 1 | F | 53 | 492 | residue_count_mismatch | ValueError('Expected 492 invariant residues, got 491') |
| 7K0J | 1 | B | 53 | 492 | residue_count_mismatch | ValueError('Expected 492 invariant residues, got 491') |
| 7K0K | 1 | B | 53 | 492 | residue_count_mismatch | ValueError('Expected 492 invariant residues, got 491') |
| 7K0M | 1 | B | 53 | 492 | residue_count_mismatch | ValueError('Expected 492 invariant residues, got 491') |
| 7K0M | 1 | H | 53 | 492 | residue_count_mismatch | ValueError('Expected 492 invariant residues, got 491') |
| 7K0N | 1 | B | 53 | 491 | residue_count_mismatch | ValueError('Expected 491 invariant residues, got 490') |
| 7K0N | 1 | H | 53 | 491 | residue_count_mismatch | ValueError('Expected 491 invariant residues, got 490') |
| 7K0O | 1 | B | 53 | 492 | residue_count_mismatch | ValueError('Expected 492 invariant residues, got 491') |
| 7K0O | 1 | H | 53 | 492 | residue_count_mismatch | ValueError('Expected 492 invariant residues, got 491') |
| 7K0P | 1 | B | 53 | 492 | residue_count_mismatch | ValueError('Expected 492 invariant residues, got 491') |
| 7K0P | 1 | F | 53 | 492 | residue_count_mismatch | ValueError('Expected 492 invariant residues, got 491') |
| 7MQ8 | 1 | WB | 46 | 245 | residue_count_mismatch | ValueError('Expected 245 invariant residues, got 193') |
| 7STX | 1 | C | 1 | 5 | residue_count_mismatch | ValueError('Expected 5 invariant residues, got 4') |
| 7T3M | 1 | C | 1 | 122 | residue_count_mismatch | ValueError('Expected 122 invariant residues, got 103') |
| 7T3M | 1 | E | 1 | 122 | residue_count_mismatch | ValueError('Expected 122 invariant residues, got 103') |
| 7T3M | 1 | F | 1 | 122 | residue_count_mismatch | ValueError('Expected 122 invariant residues, got 103') |
| 7T3M | 1 | G | 1 | 110 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 7T3M | 1 | H | 1 | 110 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 7T3M | 1 | I | 1 | 110 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 7TDW | 1 | A | 91 | 90 | residue_count_mismatch | ValueError('Expected 90 invariant residues, got 73') |
| 7UA4 | 1 | E | 2 | 107 | residue_count_mismatch | ValueError('Expected 107 invariant residues, got 103') |
| 7XOH | 1 | A | 4 | 458 | residue_count_mismatch | ValueError('Expected 458 invariant residues, got 457') |
| 7XOH | 1 | B | 4 | 458 | residue_count_mismatch | ValueError('Expected 458 invariant residues, got 457') |
| 7XOH | 1 | C | 4 | 458 | residue_count_mismatch | ValueError('Expected 458 invariant residues, got 457') |
| 7XOH | 1 | D | 4 | 458 | residue_count_mismatch | ValueError('Expected 458 invariant residues, got 457') |
| 7XYF | 1 | D | 1 | 86 | residue_count_mismatch | ValueError('Expected 86 invariant residues, got 85') |
| 7XYF | 1 | E | 1 | 106 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 7XYF | 1 | F | 1 | 105 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 7XYF | 1 | G | 1 | 95 | residue_count_mismatch | ValueError('Expected 95 invariant residues, got 82') |
| 7XYG | 1 | A | 38 | 98 | residue_count_mismatch | ValueError('Expected 98 invariant residues, got 81') |
| 7XYG | 1 | B | 38 | 98 | residue_count_mismatch | ValueError('Expected 98 invariant residues, got 80') |
| 7XYG | 1 | C | 23 | 81 | residue_count_mismatch | ValueError('Expected 81 invariant residues, got 75') |
| 7XYG | 1 | D | 18 | 86 | residue_count_mismatch | ValueError('Expected 86 invariant residues, got 74') |
| 7XYG | 1 | E | 13 | 106 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 7XYG | 1 | F | 13 | 105 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 7XYG | 1 | G | 29 | 95 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 7XYG | 1 | H | 30 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 7Y6K | 1 | C | 264 | 212 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 7YJ1 | 1 | A | 45 | 502 | residue_count_mismatch | ValueError('Expected 502 invariant residues, got 501') |
| 8ADI | 1 | F | 2 | 6 | residue_count_mismatch | ValueError('Expected 6 invariant residues, got 5') |
| 8F9X | 1 | I | 3 | 228 | residue_count_mismatch | ValueError('Expected 228 invariant residues, got 227') |
| 8FSS | 1 | B | 2 | 6 | residue_count_mismatch | ValueError('Expected 6 invariant residues, got 5') |
| 8G94 | 1 | G | 109 | 27 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 8GWG | 1 | A | 1 | 931 | residue_count_mismatch | ValueError('Expected 931 invariant residues, got 928') |
| 8GWG | 1 | C | 1 | 78 | residue_count_mismatch | ValueError('Expected 78 invariant residues, got 72') |
| 8GWG | 1 | D | 6 | 187 | residue_count_mismatch | ValueError('Expected 187 invariant residues, got 186') |
| 8GWG | 1 | I | 5 | 113 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 8I8E | 1 | A | 2 | 72 | residue_count_mismatch | ValueError('Expected 72 invariant residues, got 71') |
| 8I8E | 1 | H | 2 | 72 | residue_count_mismatch | ValueError('Expected 72 invariant residues, got 71') |
| 8I8E | 1 | I | 2 | 72 | residue_count_mismatch | ValueError('Expected 72 invariant residues, got 71') |
| 8I8E | 1 | J | 2 | 72 | residue_count_mismatch | ValueError('Expected 72 invariant residues, got 71') |
| 8I8E | 1 | K | 2 | 72 | residue_count_mismatch | ValueError('Expected 72 invariant residues, got 71') |
| 8I8E | 1 | L | 2 | 72 | residue_count_mismatch | ValueError('Expected 72 invariant residues, got 71') |
| 8IAK | 1 | B | 7 | 555 | residue_count_mismatch | ValueError('Expected 555 invariant residues, got 554') |
| 8IAK | 1 | F | 7 | 555 | residue_count_mismatch | ValueError('Expected 555 invariant residues, got 554') |
| 8PKP | 1 | S | 1 | 84 | residue_count_mismatch | ValueError('Expected 84 invariant residues, got 83') |
| 8PM2 | 1 | B | 3 | 338 | residue_count_mismatch | ValueError('Expected 338 invariant residues, got 335') |
| 8RWZ | 1 | A | 15 | 114 | residue_count_mismatch | ValueError('Expected 114 invariant residues, got 100') |
| 8RWZ | 1 | B | 1 | 160 | residue_count_mismatch | ValueError('Expected 160 invariant residues, got 146') |
| 8RWZ | 1 | E | 1 | 658 | residue_count_mismatch | ValueError('Expected 658 invariant residues, got 627') |
| 8RWZ | 1 | F | 18 | 89 | residue_count_mismatch | ValueError('Expected 89 invariant residues, got 63') |
| 8RX0 | 1 | A | 1 | 160 | residue_count_mismatch | ValueError('Expected 160 invariant residues, got 98') |
| 8RX0 | 1 | C | 1 | 97 | residue_count_mismatch | ValueError('Expected 97 invariant residues, got 91') |
| 8RX0 | 1 | H | 1 | 176 | residue_count_mismatch | ValueError('Expected 176 invariant residues, got 71') |
| 8RXA | 1 | B | 1 | 92 | invariant_not_generated | AttributeError("'NoneType' object has no attribute 'copy'") |
| 8RXO | 1 | A | 53 | 100 | residue_count_mismatch | ValueError('Expected 100 invariant residues, got 90') |
| 8SAK | 1 | E | 2 | 105 | residue_count_mismatch | ValueError('Expected 105 invariant residues, got 102') |
| 8SNI | 1 | A | 35 | 259 | residue_count_mismatch | ValueError('Expected 259 invariant residues, got 258') |
| 8VIB | 1 | C | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VIB | 1 | D | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VID | 1 | C | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 |  | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | AA | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | AD | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | AG | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | BA | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | BD | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | BG | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | C | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | CC | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | CF | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | D | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | DC | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | DF | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | EB | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | EE | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | FB | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | FE | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | GA | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | GD | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | GG | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | HA | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | HD | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | HG | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | I | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | IC | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | IF | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | J | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | JC | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | JF | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | KB | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | KE | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | LB | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | LE | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | MA | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | MD | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | MG | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | ND | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | NG | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | O | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | OC | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | OF | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | P | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | PC | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | PF | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | QB | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | QE | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | RB | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | RE | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | SA | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | SD | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | SG | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | TA | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | TD | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | TG | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | U | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | UC | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | UF | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | V | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | VC | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | VF | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | WB | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | WE | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | XB | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | XE | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | YA | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | YD | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKQ | 1 | ZA | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKQ | 1 | ZD | 45 | 93 | residue_count_mismatch | ValueError('Expected 93 invariant residues, got 85') |
| 8VKR | 1 | AA | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | AD | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | AG | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | C | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | CC | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | CF | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | EB | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | EE | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | GA | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | GD | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | GG | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | I | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | IC | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | IF | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | KB | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | KE | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | MA | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | MD | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | MG | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | O | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | OC | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | OF | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | QB | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | QE | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | SA | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | SD | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | SG | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | U | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | UC | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | UF | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | WB | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | WE | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | YA | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
| 8VKR | 1 | YD | 8 | 280 | residue_count_mismatch | ValueError('Expected 280 invariant residues, got 244') |
