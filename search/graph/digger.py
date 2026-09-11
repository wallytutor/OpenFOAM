# -*- coding: utf-8 -*-
""" Run with `uv run python digger.py` from this directory. """

import networkx as nx
from pathlib import Path
from majordome.openfoam.files import FoamCaseHandle


class FoamFvSchemesGraph:
    SELECTED_GROUPS = sorted([
        "XiFluid",
        "compressibleMultiphaseVoF",
        "compressibleVoF",
        "film",
        "fluid",
        "incompressibleDenseParticleFluid",
        "incompressibleDriftFlux",
        "incompressibleFluid",
        "incompressibleMultiphaseVoF",
        "incompressibleVoF",
        "isothermalFilm",
        "isothermalFluid",
        "movingMesh",
        "multiRegion",
        "multicomponentFluid",
        "multiphaseEuler",
        "potentialFoam",
        "shockFluid",
    ])

    __slots__ = (
        "_root",
        "_groups",
        "_max_depth",
        "_graph",
    )

    def __init__(
            self,
            tutorials_root: Path = Path("/opt/openfoam13/tutorials"),
            tutorials_groups: list[str] | None = None,
            max_depth: int = 2,
        ) -> None:
        if not tutorials_root.is_dir():
            raise NotADirectoryError(tutorials_root)

        if tutorials_groups is None:
            tutorials_groups = self.__class__.SELECTED_GROUPS

        self._root = tutorials_root
        self._groups = tutorials_groups
        self._max_depth = max_depth
        self._graph = nx.DiGraph()

    def _value_handling(self, value):
        if (key := "CrankNicolson") in value:
            value = key

        if (key := "multivariateSelection") in value:
            parts = value.replace("\n", "").split(" ")
            value = " ".join(parts[:1+parts.index(key)])

        if (key := "multivariateIndependent") in value:
            parts = value.replace("\n", "").split(" ")
            value = " ".join(parts[:1+parts.index(key)])

        return value

    @staticmethod
    def get_case_name(handle: FoamCaseHandle):
        """ Get name of tutorial with parent group. """
        parts = handle.root_dir.parts
        return "/".join(parts[1+parts.index("tutorials"):])

    def add_group(self, case_name, g, group_name, parent):
        g.add_node(group_name, node_type="group")

        for entry in parent.keys():
            value = self._value_handling(parent.get(entry))

            # Otherwise defaults become undistinguishable!
            if entry == "default":
                entry = f"default_{group_name}"

            g.add_node(entry, node_type="field")
            g.add_edge(group_name, entry)

            g.add_node(value, node_type="value")
            g.add_edge(entry, value)

            g.add_edge(value, case_name)

    def get_fv_schemes(self, handle: FoamCaseHandle):
        # TODO also test if fvSchemes* exists!
        if not (handle.root_dir / "system/fvSchemes").exists():
            # print(f"fvSchemes is missing in {handle.root_dir}")
            return

        name = self.get_case_name(handle)
        self._graph.add_node(name, node_type="case")

        def add(who, fv_schemes=handle.fv_schemes):
            self.add_group(name, self._graph, who, fv_schemes.get(who))

        add("ddtSchemes")
        add("gradSchemes")
        add("divSchemes")
        add("laplacianSchemes")
        add("snGradSchemes")

    def dig_tutorial(self, parent: Path, level: int):
        if level > self._max_depth:
            return

        for candidate in parent.iterdir():
            if not candidate.is_dir():
                continue

            case_handle = FoamCaseHandle(candidate)

            if not case_handle.is_valid:
                self.dig_tutorial(candidate, level + 1)
            else:
                self.get_fv_schemes(case_handle)

    def feed_graph(
            self,
            graph_name: str | Path | None = None,
            dump_graphml: bool = False,
            dump_gexf: bool = False,

        ):
        for group in self._groups:
            parent_dir = self._root / group

            self.dig_tutorial(parent_dir, 1)

        if not graph_name:
            return

        if isinstance(graph_name, Path):
            match graph_name.suffix:
                case ".graphml":
                    dump_graphml = True
                case ".gexf":
                    dump_gexf = True
                case _:
                    pass

            graph_name = graph_name.with_suffix("")

        if dump_graphml:
            nx.write_graphml(self._graph, f"{graph_name}.graphml")

        # Best dump for reading with Gephi
        if dump_gexf:
            nx.write_gexf(self._graph, f"{graph_name}.gexf")

    @property
    def graph(self) -> nx.DiGraph:
        return self._graph

    def subtree_from(
            self,
            start_at: str,
            selected: list[str] | None = None,
            skip_type: set[str] | None = None
        ) -> nx.DiGraph:
        skip_type = skip_type or {"case"}
        nodes_to_include = set()

        def dfs(g, start, children):
            node = g.nodes[start]

            if "node_type" in node and node["node_type"] in skip_type:
                return

            nodes_to_include.add(start)

            for child in g.successors(start):
                if children is not None and child in children:
                    dfs(g, child, None)

                if children is None:
                    dfs(g, child, None)

        dfs(self._graph, start_at, selected)

        return self._graph.subgraph(nodes_to_include).copy()


def main():
    digger =  FoamFvSchemesGraph()
    digger.feed_graph()


    groups = {
        "ddtSchemes": None,
        "gradSchemes": None,
        "divSchemes": [
            "div(phi,U)",
            "div(phi,Yi_h)",
            "div(phi,K)",
            "div(phi,k)",
            "div(phi,epsilon)",
            "div(((rho*nuEff)*dev2(T(grad(U)))))",
            # Aliases
            "energy",
            "turbulence"
        ],
        "laplacianSchemes": None,
        "snGradSchemes": None,
    }

    graphs = {}

    for name, selected in groups.items():
        graphs[name] = digger.subtree_from(name, selected=selected)
        nx.write_gexf(graphs[name], f"{name}.gexf")

    # Useful functions
    # G.successors(node)
    # G.predecessors(node)
    # nx.descendants(G, node)

    for scheme in groups["divSchemes"]:
        for kind in graphs["divSchemes"].successors(scheme):
            print(f"> {scheme}/{kind}")


if __name__ == "__main__":
    main()
