#!/usr/bin/env python

from dataclasses import dataclass
from typing import List
import json
from datetime import datetime
from collections import defaultdict
from glob import glob
import os
import subprocess
from pathlib import Path

from jinja2 import Environment, BaseLoader, Template
from bokeh.embed import components
from bokeh.models.widgets import Panel, Tabs
from bokeh.models import Div

from plotting_defaults import cyclomics_defaults, TEMPLATE_STR, nextflow_params_parser

REPORT = "report.html"


def get_template(template: str) -> Template:
    """
    get the correct template object from Jinja
    """
    # get the root directory of the project by double parenting
    # root = Path(__file__).resolve().parent.parent
    # templates_dir = root / templates_dir
    # env = Environment(loader=FileSystemLoader(templates_dir))
    # current_template = env.get_template(template)

    current_template = Environment(loader=BaseLoader()).from_string(template)
    return current_template


def human_format(num):
    if type(num) == str:
        try:
            num = float(num)
        except ValueError:
            return num

    num = float("{:.3g}".format(num))
    magnitude = 0
    while abs(num) >= 1000:
        magnitude += 1
        num /= 1000.0
    return "{}{}".format(
        "{:f}".format(num).rstrip("0").rstrip("."), ["", "K", "M", "B", "T"][magnitude]
    )


@dataclass
class SummaryCard:
    """
    icon: should be a fontawesome icon
    color: should be a bootstrap color
    """

    name: str
    icon: str
    value: str
    color: str

    def __repr__(self):
        return f"""
        <div class="col-xl-3 col-sm-6 col-12 mb-4">
            <div class="card">
                <div class="card-body">
                    <div class="d-flex justify-content-between px-md-1">
                        <div>
                            <h3 class="{self.color}">{human_format(self.value)}</h3>
                            <p class="mb-0">{self.name}</p>
                        </div>
                        <div class="align-self-center">
                            <i class="fas {self.icon} {self.color} fa-3x"></i>
                        </div>
                    </div>
                </div>
            </div>
        </div>
        """


@dataclass
class ReportTab:
    """
    name: Name of the tab in the final report
    plot: should be a Bokeh panel

    """

    name: str
    plot: str
    script: str


@dataclass
class ReportTabCollection:
    tabs: List[ReportTab]

    def __add__(self, addition):
        print(addition["name"])
        script = addition["script"].replace("\n", "")

        self.tabs.append(ReportTab(addition["name"], addition["div"], script))
        return self

    def get_scripts(self):
        return [x.script for x in self.tabs]

    def generate_tabs(self, overall_width=cyclomics_defaults.width):
        pre_tabs = []
        for tab in self.tabs:
            plot = tab.plot
            # we set the width here and it propogates to the overall width of the tab selector
            # https://github.com/bokeh/bokeh/issues/8726
            pre_tabs.append(
                Panel(child=Div(text=plot, width=overall_width), title=tab.name)
            )

        return components(Tabs(tabs=pre_tabs))


def main(args):
    html_template = get_template(TEMPLATE_STR)
    data = defaultdict(list)
    data["generation_time"] = datetime.now().strftime("%d-%b-%Y (%H:%M:%S)")
    data["additional_info"] = {}
    data["additional_info"]["git_version"] = args.version
    data["additional_info"]["nextflow_params"] = nextflow_params_parser(
        args.nextflow_params
    )

    tabs = ReportTabCollection([])

    for plot_json in glob("*.json"):
        with open(plot_json, "r") as f:
            test_json_content = json.load(f)
        # print(test_json_content)
        for k, v in test_json_content.items():
            if k == "additional_info":
                data["additional_info"].update(v)
            else:
                try:
                    tabs += v
                except (ValueError, KeyError) as e:
                    pass

    print(data["additional_info"])
    data["bokehscript"] = tabs.get_scripts()
    data["plot_items"] = tabs.generate_tabs()
    for i in [
        (
            "Sequencing reads",
            "fa-dna",
            "readsraw fastq info",
            "text-succes",
        ),  # count from input material
        (
            "Post split & QC reads",
            "fa-filter",
            "readsfiltered fastq info",
            "text-succes",
        ),
        (
            "Aligning reads",
            "fa-align-center",
            "Reference_aligned_with_backbone",
            "text-succes",
        ),  # reads aligning in initial alignment
        (
            "Consensus inserts",
            "fa-bars",
            "total_reference_mapping_reads",
            "text-succes",
        ),
        (
            "Variants found",
            "fa-map-marker-alt",
            "variants_found_non_backbone",
            "text-succes",
        ),
        ("Backbone-insert %", "fa-bullseye", "read_struc_prec_bbi", "text-succes"),
        ("Pipeline version", "fa-code-branch", "git_version", "text-succes"),
    ]:
        try:
            data["cards"].append(
                SummaryCard(i[0], i[1], data["additional_info"][i[2]], i[3])
            )
        except KeyError:
            data["cards"].append(SummaryCard(i[0], i[1], "nan", i[3]))

    # print(tabs)
    with open(REPORT, "w") as fh:
        fh.write(html_template.render(**data))


def parse_nextflow_params(params_string):
    pass


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Create hist plot from a regex for fastq and fastq.gz files."
    )

    parser.add_argument("nextflow_params", type=str)
    parser.add_argument("version", type=str)

    args = parser.parse_args()
    print(args)
    main(args)
