from abc import ABC, abstractmethod
from collections import defaultdict
from typing import Dict, List, Tuple

from .alignment import AlignmentGroup


class GroupSplitStrategy(ABC):
    @staticmethod
    @abstractmethod
    def build_alignment_splits(
        alignment_ids: List[str],
        clusters: List[Tuple[str, List[Dict], bool]],
    ) -> Dict[str, Tuple[List[str], List[Dict], bool]]:
        """
        Build mapping from cluster_id -> (matching alignment IDs, cluster_blocks, is_consecutive).

        Args:
            alignment_ids: List of alignment IDs (aln.id) to consider for grouping.
            clusters: List of tuples (cluster_id, cluster_blocks, is_consecutive).

        Returns:
            Dictionary:
                cluster_id -> (
                    matching_alignment_ids: List[str],
                    cluster_blocks: List[Dict],
                    is_consecutive: bool
                )
        """
        pass

    @abstractmethod
    def split(self, group: AlignmentGroup) -> List[AlignmentGroup]:
        """
        Splits an AlignmentGroup into one or more groups based on strategy logic.
        """
        pass


class NoSplit(GroupSplitStrategy):
    @staticmethod
    def build_alignment_splits(
        alignment_ids: List[str],
        clusters: List[Tuple[str, List[Dict], bool]],
    ) -> Dict[str, Tuple[List[str], List[Dict], bool]]:
        # No splitting — all alignments stay together
        if not clusters:
            return {}
        cluster_id, cluster_blocks, is_consecutive = clusters[0]
        return {cluster_id: (alignment_ids, cluster_blocks, is_consecutive)}

    def split(self, group: AlignmentGroup) -> List[AlignmentGroup]:
        return [group]


class RelativeIdSplit(GroupSplitStrategy):
    @staticmethod
    def build_alignment_splits(
        alignment_ids: List[str],
        clusters: List[Tuple[str, List[Dict], bool]],
    ) -> Dict[str, Tuple[List[str], List[Dict], bool]]:
        """
        Build mapping from cluster_id -> (matching alignment IDs, cluster_blocks, is_consecutive).
        """
        cluster_map = {
            cluster_id: (cluster_blocks, is_consecutive)
            for cluster_id, cluster_blocks, is_consecutive in clusters
        }

        id_to_alignments = defaultdict(list)
        for aln_id in alignment_ids:
            if aln_id in cluster_map:
                id_to_alignments[aln_id].append(aln_id)

        # Combine with cluster info
        result_map = {}
        for cluster_id, aln_ids in id_to_alignments.items():
            cluster_blocks, is_consecutive = cluster_map[cluster_id]
            result_map[cluster_id] = (aln_ids, cluster_blocks, is_consecutive)

        return result_map

    def split(
        self,
        group: AlignmentGroup,
    ) -> List[AlignmentGroup]:
        """
        Split AlignmentGroup into multiple AlignmentGroups.
        """
        # 1. Get just the IDs from the group
        alignment_ids = [aln.id for aln in group.alignments]

        # 2. Build mapping (testable in isolation)
        split_map = self.build_alignment_splits(alignment_ids, group.clusters)

        # 3. Build new groups
        result = []
        for cluster_id, (aln_ids, cluster_blocks, is_consecutive) in split_map.items():
            alignments = [aln for aln in group.alignments if aln.id in aln_ids]
            new_group = AlignmentGroup.from_existing(group, cluster_id, alignments)
            new_group.clusters = (cluster_id, cluster_blocks, is_consecutive)
            result.append(new_group)

        return result


class GroupSplitStrategyFactory:
    def get_strategy(self, structure_type: str) -> GroupSplitStrategy:
        split_keys = ["+", "chimeric", "template_switch"]
        if any(key in structure_type for key in split_keys):
            return RelativeIdSplit()

        return NoSplit()
