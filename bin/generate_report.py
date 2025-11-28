#!/usr/bin/env python

import argparse
import json
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from glob import glob
from typing import List

from jinja2 import BaseLoader, Environment, Template
from plotting_defaults import TEMPLATE_STR, human_format, nextflow_params_parser


def get_template(template: str) -> Template:
    """Load a Jinja template from string."""
    return Environment(loader=BaseLoader()).from_string(template)


@dataclass
class SummaryCard:
    name: str
    icon: str
    value: str
    color: str

    def __repr__(self):
        return f"""
        <div class="col-xl-3 col-sm-6 col-12 mb-4">
            <div class="card shadow-sm">
                <div class="card-body">
                    <div class="d-flex justify-content-between align-items-center">
                        <div>
                            <h3 class="{self.color}">{human_format(self.value)}</h3>
                            <p class="mb-0">{self.name}</p>
                        </div>
                        <i class="fas {self.icon} {self.color} fa-3x"></i>
                    </div>
                </div>
            </div>
        </div>
        """


@dataclass
class ReportTab:
    name: str
    plot: str
    script: str
    priority: int


@dataclass
class ReportTabCollection:
    tabs: List[ReportTab]

    def __add__(self, addition):
        """Allow adding tab data from JSON."""
        self.tabs.append(
            ReportTab(
                name=addition["name"],
                plot=addition.get("div", ""),
                script=addition.get("script", ""),
                priority=addition.get("priority", 0),
            )
        )
        return self

    def get_scripts(self):
        return [x.script for x in self.tabs]

    def generate_tabs(self, priority_limit=89) -> List[ReportTab]:
        """Return tabs sorted by priority."""
        return sorted(
            [x for x in self.tabs if x.priority < priority_limit],
            key=lambda t: t.priority,
        )


def main(args):
    html_template = get_template(TEMPLATE_STR)
    data = defaultdict(list)

    data["generation_time"] = datetime.now().strftime("%d-%b-%Y (%H:%M:%S)")
    data["additional_info"] = {
        "sample_name": args.sample_name,
        "git_version": args.version,
        "nextflow_params": nextflow_params_parser(args.nextflow_params),
    }

    tabs = ReportTabCollection([])

    # Load all plot .json files
    for plot_json in glob("*.json"):
        with open(plot_json, "r") as f:
            content = json.load(f)
        for k, v in content.items():
            if k == "additional_info":
                data["additional_info"].update(v)
            else:
                try:
                    tabs += v
                except (ValueError, KeyError):
                    pass

    print(data["additional_info"])
    # Tabs for rendering
    data["plotlyscript"] = tabs.get_scripts()
    data["plot_items"] = tabs.generate_tabs(priority_limit=args.priority_limit)

    # Summary cards
    data["cards"] = []
    card_definitions = [
        ("Sequencing reads", "fa-dna", "readsRaw read info", "text-dark"),
        ("Post split & QC reads", "fa-filter", "readsFiltered read info", "text-dark"),
        (
            "Aligning reads",
            "fa-align-center",
            "Reference_aligned_with_backbone",
            "text-dark",
        ),
        ("Consensus inserts", "fa-bars", "total_reference_mapping_reads", "text-dark"),
        (
            "Variants found",
            "fa-map-marker-alt",
            "variants_found_non_backbone",
            "text-dark",
        ),
        ("Pipeline version", "fa-code-branch", "git_version", "text-dark"),
    ]

    for name, icon, key, color in card_definitions:
        val = data["additional_info"].get(key, "nan")
        data["cards"].append(SummaryCard(name, icon, val, color))

    if data["additional_info"].get("exit_status") == 1:
        data["cards"].append(
            SummaryCard("Plotting error", "fa-times-circle", "Warning!", "text-warning")
        )

    report_name = (
        f"{args.sample_name}_report.html" if args.sample_name else "report.html"
    )
    with open(report_name, "w") as fh:
        fh.write(html_template.render(**data))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate Plotly report.")
    parser.add_argument("sample_name", type=str)
    parser.add_argument("nextflow_params", type=str)
    parser.add_argument("version", type=str)
    parser.add_argument("priority_limit", type=int, default=89)

    args = parser.parse_args()

    main(args)
