from collections import Counter
from dataclasses import dataclass
from typing import Dict

UNKNOWN_BARCODE = "NNNN"
BARCODE_POSITIONS = [32, 62, 79, 97]


@dataclass
class BarcodeExtrator:
    # key is Alignment, val is consensus alignment
    alignment_blocks: Dict
    offset: int = 0

    def extract_barcode(self):
        for block in self.alignment_blocks:
            if not self.is_barcodeable(block):
                result = UNKNOWN_BARCODE
                bc_array = [[]]
            else:
                result, bc_array = self.find_barcode_in_block(block, BARCODE_POSITIONS)
                break
        return result, bc_array

    def is_barcodeable(self, block, min_backbones=2, expected_length=123):
        chroms = [x.alignment_chromosome for x in list(self.alignment_blocks.keys())]

        bb_chroms = [x.startswith("BB") for x in chroms]
        bb_counts = Counter(bb_chroms)

        if bb_counts[True] >= min_backbones:
            result = True
        else:
            result = False
        return result

    def find_barcode_in_block(self, block, positions):
        bc = ""
        bc_array = []

        # Filter barcode positions that are not available due to offset
        offset_positions = [x for x in positions if x >= self.offset]

        for pos in offset_positions:
            # if we have a partial backbone we need to look earlier in the reads
            
            nucs = [
                cons.get_aligned_position(pos, self.offset)
                for cons in self.alignment_blocks.values()
            ]
            # nucs = [x for x in nucs if x != "-"]
            if len(nucs) == 0:
                # if we find nothing, break the loop and be uncertain
                bc = UNKNOWN_BARCODE
                break

            bc_array.append(nucs)
            # add most common for location
            best_nuc = Counter(nucs).most_common(1)[0][0]
            bc += best_nuc

        return bc, bc_array
