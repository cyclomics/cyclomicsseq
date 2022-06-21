from abc import ABC, abstractmethod
from collections import Counter
from typing import List

import pysam
from loguru import logger
from src.classifiers import (
    BackboneDoubleInsert,
    BackboneInsert,
    BackboneInsertGaps,
    DoubleBackbone,
    QuadtripleBackbone,
    SingleBackbone,
    SingleBackboneGaps,
    SingleDirectionalFlip,
    DirectionalFlip,
    SingleInsert,
    SingleInsertGaps,
    ShortContiguous,
    ShortGaps,
    TripleBackbone,
    TripleInsert,
    TwoInsertsChromosomal,
    TwoInsertsChromosomalGaps,
    TwoInsertsStartLocations,
    Singleton,
    SingletonBackbone,
    ComplexConcatemer,
)


class BaseClassifier(ABC):
    def __init__(self):
        self.result = {}
        self.result_details = {}

    @abstractmethod
    def classify(self, group) -> str:
        return "Read_type"

    def create_piechart(self):
        pass

    def show_counts(self):
        return Counter(self.result.values())

    def get_metadata(self, classifier_name, group, consensus):
        classifier = [x for x in self.classifiers if x.name == classifier_name][0]
        metadata = classifier.create_metadata(group, consensus)
        return metadata


class ClassifierGroup(BaseClassifier):
    def __init__(self):
        logger.debug(f"Classifier group enabled: {self.__class__.__name__}")
        # define all classifiers to use, priority: lower is more important
        self.classifiers = [
            ShortGaps(priority=99),
            ShortContiguous(priority=98),
            Singleton(priority=95),
            SingletonBackbone(priority=94),
            ComplexConcatemer(priority=93),
            QuadtripleBackbone(priority=24),
            SingleInsertGaps(priority=22),
            SingleInsert(priority=18),
            SingleBackboneGaps(priority=17),
            SingleBackbone(priority=16),
            SingleDirectionalFlip(priority=15),
            DirectionalFlip(priority=14),
            TripleInsert(priority=12),
            TripleBackbone(priority=11),
            BackboneDoubleInsert(priority=10),
            TwoInsertsChromosomalGaps(priority=8),
            TwoInsertsStartLocations(priority=7),
            TwoInsertsChromosomal(priority=6),
            DoubleBackbone(priority=5),
            BackboneInsertGaps(priority=2),
            BackboneInsert(priority=1),
        ]
        self.classifiers.sort(key=lambda x: x.priority)
        logger.debug("done sorting")
        super().__init__()

    def classify(self, name, group) -> List[str]:
        logger.debug(f"Classifier group called: {self.__class__.__name__}")
        if not group:
            raise IOError

        # if we dont have any data, the result is simple
        if len(group) == 0:
            result = "Filtered"
            return [result, ""]

        outcomes = [x.classify(group) for x in self.classifiers]
        # filter on succes state
        filtered_outcomes = [x for x in outcomes if x[2]]

        # results are ordered by priority, so first succes is the best
        if filtered_outcomes:
            best_outcome = filtered_outcomes[0]
            self.result[name] = best_outcome[0]
        else:
            best_outcome = (None, None)
            self.result[name] = "Unkown"

        # Update the read specifics for detailed output
        read_specifics = {}
        for classifier_outcome in outcomes:
            read_specifics[classifier_outcome[0]] = classifier_outcome[-1]
            read_specifics[classifier_outcome[0]]["passed"] = classifier_outcome[2]

            self.result_details[name] = read_specifics

        # leave the details here
        return best_outcome[:2]


class SimpleClassifierGroup(BaseClassifier):
    def __init__(self):
        logger.debug(f"Classifier group enabled: {self.__class__.__name__}")
        self.classifiers = [ShortGaps(priority=99), BackboneInsertGaps(priority=2)]
