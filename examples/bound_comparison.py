import itertools as it
import sys

from dataclasses import dataclass
from typing import (
    Callable,
    List,
    Tuple,
)

# Sage
import sage.all

from sage.graphs.graph import Graph  # type: ignore
from sage.graphs import graph_generators
graph_generators = graph_generators.GraphGenerators()

# Local
from association_schemes.bounds import hoffman_bound, lp_bound

def actual_bound(graph: Graph) -> int:
    return len(graph.independent_set())

bounds = {
    "hoffman": hoffman_bound,
    "lp": lp_bound,
    "actual": actual_bound,
}

def _make_graph(
    s: str,
    G: Callable[..., Graph],
) -> Callable[[Tuple], Tuple[str, Graph]]:
    def func(vals: Tuple) -> Tuple[str, Graph]:
        return (
            s + str(vals),
            G(*vals)
        )
    return func

graphs = dict(
    list(map(_make_graph("H", graph_generators.HammingGraph), [
        (4, 2), (4, 3),
        (5, 2), (5, 3),
    ]))
    + list(map(_make_graph("J", graph_generators.JohnsonGraph), [
        (6, 2), (6, 3),
        (7, 2), (7, 3),
    ]))
    + list(map(_make_graph("G", graph_generators.GrassmannGraph), [
        (2, 4, 2), (2, 5, 2),
        (3, 4, 2),
    ]))
)

experiments = [
    (bound, graph)
    for bound, graph in it.product(bounds.keys(), graphs.keys())
    if bounds != "actual" or not graph.startswith("G")
]

def run_experiments(experiments):
    for bound, graph in experiments:
        print(bound, graph, file=sys.stderr)
        yield (
            bound,
            graph,
            bounds[bound](graphs[graph]),
        )

if __name__ == "__main__":
    print("bound,graph,value")
    for bound, graph, value in run_experiments(experiments):
        print(f"{bound},\"{graph}\",{value}")
