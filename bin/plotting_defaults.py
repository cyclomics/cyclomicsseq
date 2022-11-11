from dataclasses import dataclass
import re


@dataclass
class plot_defaults:
    width: int = 1100
    plot_title_size: str = "18pt"
    plot_axis_text_size: str = "16pt"
    plot_label_size: str = "12pt"


def findnth_occurance(input, search, n):
    """
    Finding nth occurrence of search in input, zero indexed.
    returns -1 if no nth occurance exists

    abcabcabc, c, 0 -> 2
    """
    #
    inilist = [m.start() for m in re.finditer(search, input)]
    # correct length for 0 indexing iso 1 indexing
    if len(inilist) > n:
        return inilist[n]
    else:
        return -1


def nextflow_params_parser(params_str):
    """
    recursive function to turn the params string from nextflow into a dictionary.

    eg turn:
    [a:b, c:[d:e, f:g, h:[i:j], h2:[i2:[nested:3]]], k:l,m:n]
    into a dict like:
    {
        'a': 'b',
        'c': {
            'd': 'e',
            'f': 'g',
            'h': {'i': 'j'},
            'h2': {
                'i2': {
                    'nested': '3'
                }
            }
        },
        'k': 'l',
        'm': 'n'
    }
    """
    # initially drop the first and last brackets
    result = dict()

    # drop containment sqr brackets, remove whitespace before
    params_str = params_str.rstrip()[1:-1]
    while "," in params_str or "[" in params_str:
        # get key by finding secondary separator
        next_sep = params_str.find(":")
        key = params_str[:next_sep]

        # if we find a nested list we rerun this function on the nested list
        if params_str[next_sep + 1] == "[":
            closer_count = 0
            closer_idx = findnth_occurance(params_str, "]", closer_count)

            substring = params_str[params_str.find("[") : closer_idx + 1]
            # count openers in this substring:
            opens = substring.count("[")
            closes = substring.count("]")
            while opens != closes:
                closer_count += 1
                closer_idx = findnth_occurance(params_str, "]", closer_count)
                substring = params_str[params_str.find("[") : closer_idx + 1]
                opens = substring.count("[")
                closes = substring.count("]")

            value = nextflow_params_parser(substring)
            item_end = closer_idx + 1

        # otherwise we find the next item in the list, and work our way back to find the last comma
        else:
            next_next_sep = params_str[next_sep + 1 :].find(":")
            item_end = params_str[: next_sep + next_next_sep].rfind(",")
            value = params_str[next_sep + 1 : item_end]

        # remove the observed key and value, and save them
        params_str = params_str[item_end + 1 :].lstrip()
        result.update({key: value})

    # get the last item if available, not true if nested list with one items
    if params_str:
        next_sep = params_str.find(":")
        key = params_str[:next_sep]
        value = params_str[next_sep + 1 :]
        result.update({key: value})

    return result


cyclomics_defaults = plot_defaults()

# We cannot guarantee the presence of a file in the nextflow container
TEMPLATE_STR = """
<!doctype html>
<html lang="en">
   <head>
      <!-- Required meta tags -->
      <meta charset="utf-8">
      <meta name="viewport" content="width=device-width, initial-scale=1, shrink-to-fit=no">
      <!-- Bootstrap -->
      <link rel="stylesheet" href="https://stackpath.bootstrapcdn.com/bootstrap/4.3.1/css/bootstrap.min.css" integrity="sha384-ggOyR0iXCbMQv3Xipma34MD+dH/1fQ784/j6cY/iJTQUOhcWr7x9JvoRxT2MZw1T" crossorigin="anonymous">
      <!-- Fontawesome -->
      <script src="https://kit.fontawesome.com/3282f62a28.js" crossorigin="anonymous"></script>
      <title>
         Cyclomics Report
      </title>
   </head>
   <body>
      <div class="container body">
         <div class="main_container">
            <div class="container">
               <nav class="navbar navbar-light bg-light my-4 border-bottom border-dark">
                  <a class="navbar-brand" href="#">
                     <!-- logo and Name -->
                     <img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAJYAAACWCAYAAAA8AXHiAAAACXBIWXMAAC4jAAAuIwF4pT92AAAGvmlUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPD94cGFja2V0IGJlZ2luPSLvu78iIGlkPSJXNU0wTXBDZWhpSHpyZVN6TlRjemtjOWQiPz4gPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iQWRvYmUgWE1QIENvcmUgNS42LWMxNDAgNzkuMTYwNDUxLCAyMDE3LzA1LzA2LTAxOjA4OjIxICAgICAgICAiPiA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPiA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIiB4bWxuczp4bXA9Imh0dHA6Ly9ucy5hZG9iZS5jb20veGFwLzEuMC8iIHhtbG5zOnhtcE1NPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvbW0vIiB4bWxuczpzdEV2dD0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL3NUeXBlL1Jlc291cmNlRXZlbnQjIiB4bWxuczpkYz0iaHR0cDovL3B1cmwub3JnL2RjL2VsZW1lbnRzLzEuMS8iIHhtbG5zOnBob3Rvc2hvcD0iaHR0cDovL25zLmFkb2JlLmNvbS9waG90b3Nob3AvMS4wLyIgeG1wOkNyZWF0b3JUb29sPSJBZG9iZSBQaG90b3Nob3AgQ0MgKE1hY2ludG9zaCkiIHhtcDpDcmVhdGVEYXRlPSIyMDE5LTAxLTMwVDEwOjU0OjQyKzAxOjAwIiB4bXA6TWV0YWRhdGFEYXRlPSIyMDE5LTAxLTMwVDEzOjQ4OjE3KzAxOjAwIiB4bXA6TW9kaWZ5RGF0ZT0iMjAxOS0wMS0zMFQxMzo0ODoxNyswMTowMCIgeG1wTU06SW5zdGFuY2VJRD0ieG1wLmlpZDpkN2Y5MjQ2ZS1jMzEzLTRjOGYtOWQwYi0yZmU0MGRkZWMwYWEiIHhtcE1NOkRvY3VtZW50SUQ9ImFkb2JlOmRvY2lkOnBob3Rvc2hvcDo4YjgyNTE4Zi01ZTFhLTFiNDUtOTU3ZC01ZjFiMWE1ZGNhOWMiIHhtcE1NOk9yaWdpbmFsRG9jdW1lbnRJRD0ieG1wLmRpZDpkMTVjODRjZC0yNmM5LTRiNDctYjM0YS00ZTRjN2EwMmY4OTMiIGRjOmZvcm1hdD0iaW1hZ2UvcG5nIiBwaG90b3Nob3A6Q29sb3JNb2RlPSIzIiBwaG90b3Nob3A6SUNDUHJvZmlsZT0ic1JHQiBJRUM2MTk2Ni0yLjEiPiA8eG1wTU06SGlzdG9yeT4gPHJkZjpTZXE+IDxyZGY6bGkgc3RFdnQ6YWN0aW9uPSJjcmVhdGVkIiBzdEV2dDppbnN0YW5jZUlEPSJ4bXAuaWlkOmQxNWM4NGNkLTI2YzktNGI0Ny1iMzRhLTRlNGM3YTAyZjg5MyIgc3RFdnQ6d2hlbj0iMjAxOS0wMS0zMFQxMDo1NDo0MiswMTowMCIgc3RFdnQ6c29mdHdhcmVBZ2VudD0iQWRvYmUgUGhvdG9zaG9wIENDIChNYWNpbnRvc2gpIi8+IDxyZGY6bGkgc3RFdnQ6YWN0aW9uPSJzYXZlZCIgc3RFdnQ6aW5zdGFuY2VJRD0ieG1wLmlpZDo3YjBmNTdkYS1iNGM4LTQ2MTEtYjBlYi1lM2Y4MTc4ZWI2MTAiIHN0RXZ0OndoZW49IjIwMTktMDEtMzBUMTA6NTQ6NDIrMDE6MDAiIHN0RXZ0OnNvZnR3YXJlQWdlbnQ9IkFkb2JlIFBob3Rvc2hvcCBDQyAoTWFjaW50b3NoKSIgc3RFdnQ6Y2hhbmdlZD0iLyIvPiA8cmRmOmxpIHN0RXZ0OmFjdGlvbj0ic2F2ZWQiIHN0RXZ0Omluc3RhbmNlSUQ9InhtcC5paWQ6ZDdmOTI0NmUtYzMxMy00YzhmLTlkMGItMmZlNDBkZGVjMGFhIiBzdEV2dDp3aGVuPSIyMDE5LTAxLTMwVDEzOjQ4OjE3KzAxOjAwIiBzdEV2dDpzb2Z0d2FyZUFnZW50PSJBZG9iZSBQaG90b3Nob3AgQ0MgKE1hY2ludG9zaCkiIHN0RXZ0OmNoYW5nZWQ9Ii8iLz4gPC9yZGY6U2VxPiA8L3htcE1NOkhpc3Rvcnk+IDwvcmRmOkRlc2NyaXB0aW9uPiA8L3JkZjpSREY+IDwveDp4bXBtZXRhPiA8P3hwYWNrZXQgZW5kPSJyIj8+MtjVngAASqdJREFU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"
                        width="30" height="30" class="d-inline-block align-top" alt="Cyclomics logo">
                     CyclomicsSeq Amplicon RUO kit Report
                  </a>
               </nav>
            </div>
            <div class="row">
               <div class="col-12 mt-3 mb-1">
                  <h5 class="text-uppercase">Quick information</h5>
                  <!-- <p>Statistics on minimal cards.</p> -->
               </div>
            </div>
            <div class="row">
               {% for card in cards %}
               {{ card }}
               {% endfor %}
            </div>
            <div class="col-12 mt-3 mb-1">
               <h3>Inputs</h3>
               <div id="accordion">
                  <div class="card">
                     <div class="card-header" id="headingOne">
                        <h5 class="mb-0">
                           <button class="btn btn-link collapsed" data-toggle="collapse" data-target="#collapseOne" aria-expanded="true" aria-controls="collapseOne">
                           Details and Inputs
                           </button>
                        </h5>
                     </div>
                     <div id="collapseOne" class="collapse" aria-labelledby="headingOne" data-parent="#accordion">
                        <div class="card-body">
                           <ul class="list-group"> 
                              <li class="list-group-item">Read directory: {{ additional_info.nextflow_params.input_read_dir }}</li>
                              <li class="list-group-item">Read pattern: {{ additional_info.nextflow_params.read_pattern }}</li>
                              <li class="list-group-item">Reference: {{ additional_info.nextflow_params.reference }}</li>
                              <li class="list-group-item">Backbone: {{ additional_info.nextflow_params.backbone }}</li>
                              <li class="list-group-item">Regions: {{  additional_info.nextflow_params.region_file }}</li>
                              <li class="list-group-item">Original output: {{  additional_info.nextflow_params.output_dir }}</li>
                              &nbsp
                              <li class="list-group-item">QC: {{  additional_info.nextflow_params.qc }}</li>
                              <li class="list-group-item">Consensus generation: {{  additional_info.nextflow_params.consensus_calling }}</li>
                              <li class="list-group-item">Alignment strategy: {{  additional_info.nextflow_params.alignment }}</li>
                              <li class="list-group-item">Variant calling: {{  additional_info.nextflow_params.variant_calling }}</li>
                              &nbsp
                              <li class="list-group-item">Runtime environment: {{  additional_info.nextflow_params.profile_selected }}</li>
                              <li class="list-group-item">generation date: {{ generation_time }}</li>
                              <li class="list-group-item">Software version: {{ additional_info.git_version }}</li>
                           </ul>
                        </div>
                     </div>
                  </div>
               </div>
            </div>
            <div class="row" id="plots">
               <div class="col-12 mt-5 mb-1">
                  <h3>Plots and Graphs</h3>
                  {% for plot_item in plot_items %}
                  {{ plot_item }}
                  {% endfor %}
               </div>
            </div>
            <!-- need to go below JS for layout reasons -->
            {% for i in bokehscript %}
            {{ i }}
            {% endfor %}
         </div>
      </div>
      <script src="https://code.jquery.com/jquery-3.2.1.slim.min.js" integrity="sha384-KJ3o2DKtIkvYIK3UENzmM7KCkRr/rE9/Qpg6aAZGJwFDMVNA/GpGFF93hXpG5KkN" crossorigin="anonymous"></script>
      <script src="https://cdnjs.cloudflare.com/ajax/libs/popper.js/1.12.9/umd/popper.min.js" integrity="sha384-ApNbgh9B+Y1QKtv3Rn7W3mgPxhU9K/ScQsAP7hUibX39j7fakFPskvXusvfa0b4Q" crossorigin="anonymous"></script>
      <script src="https://maxcdn.bootstrapcdn.com/bootstrap/4.0.0/js/bootstrap.min.js" integrity="sha384-JZR6Spejh4U02d8jOt6vLEHfe/JQGiRRSQQxSfFWpi1MquVdAyjUar5+76PVCmYl" crossorigin="anonymous"></script>
      <!-- TODO: ADD bootstrap van fontawesome js-->
      <script src="https://cdn.bokeh.org/bokeh/release/bokeh-2.4.0.min.js"
         crossorigin="anonymous"></script>
      <script src="https://cdn.bokeh.org/bokeh/release/bokeh-widgets-2.4.0.min.js"
         crossorigin="anonymous"></script>
      <script src="https://cdn.bokeh.org/bokeh/release/bokeh-tables-2.4.0.min.js"
         crossorigin="anonymous"></script>
      <script src="https://cdn.bokeh.org/bokeh/release/bokeh-gl-2.4.0.min.js"
         crossorigin="anonymous"></script>
      <script src="https://cdn.bokeh.org/bokeh/release/bokeh-mathjax-2.4.0.min.js"
         crossorigin="anonymous"></script>
      <script>
         function sleep(ms) {
           return new Promise(resolve => setTimeout(resolve, ms));
         }
         
         document.getElementById('plots').style.display = 'none';
         sleep(200)
         document.getElementById('plots').style.display = 'block';
         
      </script>
   </body>
</html>
"""
