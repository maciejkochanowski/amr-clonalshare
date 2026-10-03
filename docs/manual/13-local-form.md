# 13 The local form, and the Windows folder

The package has one engine and two ways to drive it. The command line
(`amr-clonalshare --config ...`) is the one the manual describes, and the one
the records of this release were written with. The local form drives the same
engine from the browser, for a laboratory that works without a command line,
and writes the same configuration file, the same record and the same reports.

## Start the form

```bash
amr-clonalshare-gui
```

The command starts a small web server on the loopback address of this computer
(127.0.0.1, a free port), prints the address, and opens the browser on it. The
address carries a token of the session; a page opened without it is refused,
so no other programme on the computer, and no web page, can drive a run. The
server answers this computer only and sends nothing anywhere. `--no-browser`
prints the address without opening the browser, `--port` fixes the port,
`--threads` limits the processor threads of a run, and `--root` names the
working folder (default: `AMR-ClonalShare` in the user folder;
`AMR_CLONALSHARE_HOME` in the environment sets the same). Close the window,
or press Ctrl+C, to stop.

## One run

The first page takes the tables: the metadata table, the call table, the
MIC table, as files or as cells pasted from a spreadsheet under "Paste a
table instead of choosing a file". The metadata table may be left out when
the call or the MIC table carries the lineage column itself: one sheet with
the identifier, the lineage and one column per agent is then a complete
input. The second page shows the columns of each table in lists, with the
guess of the same rules `--init` uses already chosen: the isolate identifier,
the lineage, the intake, the strata of the calls, the sampling unit, the
contrast; for the calls the antimicrobial and the call columns of a long
table or the agent columns of a wide one, the kind of call and the duplicate
policy; for the MIC table the antimicrobial and measurement columns of a long
table or the agent columns of a wide one, the censoring sign, the unit, the
laboratory and the covariates. Agent codes (`CIP`, `CIP_NM`) are shown with
the names they are read as. Under "More" of the MIC section the recorded
concentrations of every agent are listed with the share of readings on the
lowest or highest of them; "Use the observed range" writes them into the
tested-concentrations field, a panel preset fills the field from the
published range of an EU harmonised panel, and a warning asks for the range
where many readings sit on the two ends while none is declared. The budgets,
the seed and the reading of the end wells follow. Every choice is a guess to
be read, not a decision the form takes.

"Check the input only" writes the input check and stops, which is the right
first run on a new table. "Run the analysis" starts the run in the
background; the page of the run refreshes itself, shows the progress lines
the command line prints, and when the run ends lists the results: the report,
the results table, the record, the input check, and everything as one zip
archive. A table of the lineage shares, read from the record, is shown above
the files, with the conclusion for every agent in its last column (established,
not established, whole range, not estimable) and the sentence behind each
conclusion folded below the table, together with a
paragraph for a laboratory report that one button copies. "Open the folder of
this run" opens it in the file manager; "Repeat with new tables" opens the
New run page with the column choices, the tested concentrations, the budgets
and the seed of the run proposed for the next tables. "Run the analysis"
reads the tables once before the run starts, so a column named wrongly comes
back to the Columns page with its reason rather than ending the run on its
page. The "Runs" page lists every run of the working folder and sets a run
aside on request: its folder moves under `_removed` beside the others, and
nothing is deleted. Every page ends with a button that copies its address,
and with a link to the manual when the folder ships one. Closing the console
window stops the server; Ctrl+C during a run asks once more before
stopping it.

Every run has a folder of its own under the working folder, named by the time
it started and its name: the uploaded tables under `data/`, the configuration
the form wrote as `config.yaml`, the results under `results/`, and a log. The
run is therefore repeatable from the command line, on the same computer or
another, with `amr-clonalshare --config <folder>/config.yaml --results-dir
<new folder> --seed <seed>`; the page of the run prints that command with
the seed of the run. The "Runs" page lists
the runs of the working folder, those of earlier sessions included.

A few settings stay with the configuration file, because they are rarely
changed: the units per agent (`mic_units`, rarely needed when a unit column
exists), the folds and fold assignments of the e-values (`evidence.folds`,
`evidence.repeats`; the form's folds set those of the share), the level and
gates of the decomposition (`surveillance.q_fdr`,
`surveillance.min_shared_support`, `surveillance.label_alpha`) and the
`enabled` switch of each section. [The input manual](02-input.md) describes
the units and [the run page](04-run.md) the other settings.

## The examples, and a configuration file of your own

The "Examples" page lists the example collections of the repository: every
configuration file one folder below `examples/` (the fictional 96-isolate
collection of `examples/workflows`, which runs in seconds; the *E. coli*
collection, which takes seconds; the *S. suis* and the *Salmonella*
collections, which take a few minutes), each with the opening comment of
its file, and the matched tables of `examples/matched_lineages` and
`examples/ecoli_swine` for the comparison. One click runs a collection
or checks its input: the tables and the configuration are copied into a
folder of the run, so the run is the one the command line makes with
`amr-clonalshare --config examples/<collection>/<file>.yaml`, and it is
repeated from the folder like any other. A source checkout finds its
`examples` folder by itself; an installed wheel carries none, so
`--examples <folder>` names one (the Windows folder passes its own).

The same page runs a configuration file already on this computer, as
`amr-clonalshare --config` does: give the full path of the file, and the run
reads the tables where the file names them, relative to the file, and gets a
folder of its own for its configuration and its results.

## Two lineage definitions

The "Compare" page runs `amr-clonalshare-compare` on one table with one row
per isolate: the identifier, a 0/1 outcome and the two lineage columns. Its
results are the same four files the command line writes, and a second comparison of the same table goes to a folder of its own beside the first.

## The Windows folder

For a computer without Python, every release carries
`AMR-ClonalShare-1.0.0-windows.zip`, attached to the release on GitHub. It
holds the embeddable distribution of Python from python.org, the runtime
part of the pinned stack of `requirements-lock.txt` as Windows wheels, the
package, the `examples` folder of the repository, the manual as HTML pages
(`manual/`, served by the form under its "Manual" link), a two-page
`QUICK_START.pdf`, and two files to start from: `AMR-ClonalShare.bat` for
the form and `amr-clonalshare.cmd` for the command line
(`amr-clonalshare.cmd example ecoli --results-dir <folder>`). Unpack the zip
anywhere and double-click the first: a window opens and stays open while you
work, and the browser opens on the form; the window prints the address to
paste if the browser does not open. Windows may warn before the first start
that the folder is not from a known publisher (SmartScreen): choose "More
info", then "Run anyway"; the folder is not signed, and its SHA-256 is beside
the zip. Nothing is installed, no administrator rights are needed, and the
folder is removed by deleting it.

The folder is built by `scripts/build_windows_bundle.py` from the same wheels
as the release, on any operating system, and checked on a Windows runner at
every push (`.github/workflows/windows.yml`): the shipped *S. suis*
configuration run with the folder's interpreter must reproduce the shipped
record, and the tests of the form must pass with it. The interpreter of the
folder is Python 3.13, not the 3.11 of the records; the numbers agree to the
round-off the reproduction guide describes, and the check requires it.
