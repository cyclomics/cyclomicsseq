import logging
import time
from collections import Counter
from math import pi

import pandas as pd
from bokeh.io import output_file, save
from bokeh.plotting import figure
from bokeh.transform import cumsum


from src.classifiers import BaseClassifier


# Create a dict of colors for each posible classifier
CLASS_COLORS = {}
for v in BaseClassifier.__subclasses__():
    initialized = v(0)
    CLASS_COLORS[initialized.name] = initialized.color
    CLASS_COLORS["Unkown"] = "#D3D3D3"
    CLASS_COLORS["Filtered"] = "#A9A9A9"


def create_piechart(classes: dict, title: str = "Read classification") -> None:
    classes = Counter(classes.values())
    data = (
        pd.Series(classes)
        .reset_index(name="value")
        .rename(columns={"index": "read_type"})
    )
    data["angle"] = data["value"] / data["value"].sum() * 2 * pi
    data["color"] = [CLASS_COLORS[x] for x in data["read_type"]]

    p = figure(
        height=500,
        width=800,
        title=title,
        toolbar_location=None,
        tools="hover",
        tooltips="@read_type: @value",
        x_range=(-0.5, 1.0),
    )

    p.wedge(
        x=0,
        y=1,
        radius=0.4,
        start_angle=cumsum("angle", include_zero=True),
        end_angle=cumsum("angle"),
        line_color="white",
        fill_color="color",
        legend_field="read_type",
        source=data,
    )

    p.axis.axis_label = None
    p.axis.visible = False
    p.grid.grid_line_color = None

    timestamp = int(time.time())

    piechart_filename = f"plots/{timestamp}_piechart.html"

    output_file(filename=piechart_filename, title="piechart of identified read classes")
    save(p)
    logging.info(f"Created plot at {piechart_filename}")
