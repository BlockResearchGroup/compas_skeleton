from math import radians

from compas.datastructures import Datastructure
from compas.datastructures import Graph
from compas.datastructures import Mesh
from compas.datastructures.graph.duality import node_sort_neighbors
from compas.geometry import Line
from compas.geometry import Polyline
from compas.geometry import intersection_line_line
from compas.geometry import offset_polyline
from compas.itertools import flatten
from compas.itertools import pairwise
from compas.itertools import window


class SkeletonGraph(Graph):
    """The base graph of a skeleton pattern.

    Attributes
    ----------
    node_width : float
        The width of the coarse mesh at the graph nodes.
    leaf_angle : float
        The angle formed by the boundary edges of the coarse mesh
        connected to the vertices corresponding to leaf nodes of the graph.
    loops : list[list[int]]
        A list of loops formed in between leaf nodes.
        Each loop is an ordered list of nodes from one leaf to another.
    polylines : list[list[point]]
        The offset polylines corresponding to the loops.
        The offsets are used to create the topology and initial geometry of the coarse mesh.

    """

    def __init__(self):
        super().__init__()
        self.update_default_node_attributes(neighbors=None)
        self._node_width = 0
        self._leaf_angle = 0
        self._loops: list[list[int]] = []
        self._polylines: list[Polyline] = []

    @property
    def node_width(self) -> float:
        if not self._node_width:
            number_of_edges = self.number_of_edges()
            length = sum(self.edge_length(e) for e in self.edges())
            self._node_width = 0.3 * length / number_of_edges
        return self._node_width

    @node_width.setter
    def node_width(self, value: float) -> None:
        self._node_width = value

    @property
    def leaf_angle(self) -> float:
        return self._leaf_angle

    @leaf_angle.setter
    def leaf_angle(self, value: float) -> None:
        self._leaf_angle = value

    @property
    def loops(self) -> list[list[int]]:
        if not self._loops:
            self._loops = []
            for leaf in self.leaves():
                nbr = self.neighbors(leaf)[0]
                self._loops.append(self.find_edge_loop([leaf, nbr]))
        return self._loops

    @property
    def polylines(self) -> list[Polyline]:
        if not self._polylines:
            self._polylines = []
            for loop in self.loops:
                points = [self.node_coordinates(node) for node in loop]
                offset = offset_polyline(points, distance=self.node_width)
                self._polylines.append(Polyline(offset))
        return self._polylines

    @classmethod
    def from_lines(cls, lines, node_width=None, leaf_angle=None, precision=None):
        graph = super().from_lines(lines, precision)
        if node_width:
            graph.node_width = node_width
        if leaf_angle:
            graph.leaf_angle = leaf_angle
        graph.order_cycles()
        return graph

    def order_cycles(self, ccw: bool = True) -> None:
        node_xyz = {node: self.node_point(node) for node in self.nodes()}
        for node in self.nodes():
            nbrs = self.neighbors(node)
            if not self.is_leaf(node):
                nbrs = node_sort_neighbors(node, nbrs, node_xyz, ccw=ccw)
            self.node_attribute(node, name="neighbors", value=nbrs)

    def find_edge_loop(self, edge: tuple[int, int]) -> list[int]:
        u, v = edge
        loop = [u]
        while True:
            loop.append(v)
            if self.is_leaf(v):
                break
            nbrs: list[int] = self.node_attribute(v, name="neighbors")[::-1]
            nbr = nbrs[nbrs.index(u) - 1]
            u, v = v, nbr
        return loop


class SkeletonMesh(Mesh):
    """The coarse mesh of a skeleton pattern.

    Attributes
    ----------
    loop_vertices : dict[tuple[int], tuple[int]]
        A dictionary mapping loops in the skeleton graph to corresponding loops
        of boundary vertices in the coarse mesh.

    """

    def __init__(
        self,
        default_vertex_attributes=None,
        default_edge_attributes=None,
        default_face_attributes=None,
        name=None,
        **kwargs,
    ):
        super().__init__(
            default_vertex_attributes,
            default_edge_attributes,
            default_face_attributes,
            name,
            **kwargs,
        )
        self.update_default_vertex_attributes(node=None, triplet=None)
        self.update_default_edge_attributes(crease=0)
        self._loop_vertices = {}

    @property
    def loop_vertices(self) -> dict[tuple, tuple]:
        return self._loop_vertices

    @classmethod
    def from_graph(cls, graph: SkeletonGraph) -> "SkeletonMesh":
        mesh = cls()

        for node in graph.nodes():
            x, y, z = graph.node_coordinates(node)
            vertex = mesh.add_vertex(key=node, x=x, y=y, z=z, node=node)

        loop_vertices = {}

        for loop, polyline in zip(graph.loops, graph.polylines):
            triplets = list(window(loop[:1] + loop + loop[-1:], 3))
            vertices = []
            for (u, v, w), (x, y, z) in zip(triplets, polyline):
                vertex = mesh.add_vertex(x=x, y=y, z=z, triplet=(u, v, w))
                vertices.append(vertex)
            loop_vertices[tuple(loop)] = vertices

        for loop in loop_vertices:
            vertices = loop_vertices[loop]
            for (u, v), (uu, vv) in zip(pairwise(loop), pairwise(vertices)):
                mesh.add_face([uu, u, v, vv])

        if graph.leaf_angle != 0:
            angle = radians(graph.leaf_angle)

            for loop, vertices in loop_vertices.items():
                u = loop[0]
                uu, vv = vertices[0:2]
                a, b, d = mesh.vertices_attributes(names="xyz", keys=[u, uu, vv])
                ab = Line(a, b)
                bd = Line(b, d)
                ab.rotate(angle=-angle, axis=[0, 0, 1], point=a)
                x = intersection_line_line(ab, bd)[0]
                if x:
                    mesh.vertex_attributes(uu, names="xyz", values=x)

                u = loop[-1]
                vv, uu = vertices[-2:]
                a, b, d = mesh.vertices_attributes(names="xyz", keys=[u, uu, vv])
                ab = Line(a, b)
                bd = Line(b, d)
                ab.rotate(angle=angle, axis=[0, 0, 1], point=a)
                x = intersection_line_line(ab, bd)[0]
                if x:
                    mesh.vertex_attributes(uu, names="xyz", values=x)

        mesh._loop_vertices = loop_vertices
        return mesh


class Skeleton(Datastructure):
    """Compotise data structure combining the skeleton graph and coarse mesh to produce skeleton patterns.

    Parameters
    ----------
    lines : list[:class:`Line`]
        The input lines.

    Attributes
    ----------
    node_width : float
        Width of the coarse mesh at the nodes of the graph.
    leaf_angle : float
        The angle formed by the boundary edges of the coarse mesh connected to the graph leaves.
    density : int
        Number of subd iterations used to generate the final skeleton pattern.
    graph : :class:`SkeletonGraph`
        The graph used to manage the overall topology of the pattern.
    mesh : :class:`SkeletonMesh`
        The coarse mesh from which the final pattern is generated using (creased) subdivision.
    pattern : :class:`Mesh`
        The final skeleton mesh pattern.

    """

    def __init__(self, lines: list[Line]):
        super().__init__()
        self._density = 1
        self._graph = SkeletonGraph.from_lines(lines)
        self._mesh = None
        self._mesh_is_up_to_date = False
        self._pattern = None

    @property
    def node_width(self) -> float:
        return self.graph.node_width

    @node_width.setter
    def node_width(self, value: float) -> None:
        self.graph.node_width = value
        self._mesh_is_up_to_date = False
        self._pattern = None

    @property
    def leaf_angle(self) -> float:
        return self.graph.leaf_angle

    @leaf_angle.setter
    def leaf_angle(self, value: float) -> None:
        self.graph.leaf_angle = value
        self._mesh_is_up_to_date = False
        self._pattern = None

    @property
    def density(self) -> int:
        return self._density

    @density.setter
    def density(self, value: int) -> None:
        self._density = value
        self._pattern = None

    @property
    def graph(self) -> SkeletonGraph:
        return self._graph

    @property
    def mesh(self) -> SkeletonMesh:
        if self._mesh is None:
            self.init_mesh()
            self._mesh_is_up_to_date = True
        if not self._mesh_is_up_to_date:
            self.update_mesh()
        return self._mesh

    @property
    def pattern(self) -> Mesh:
        if self._pattern is None:
            self.compute_pattern()
        return self._pattern

    def init_mesh(self) -> None:
        self._mesh = SkeletonMesh.from_graph(self.graph)

    def update_mesh(self) -> None:
        for vertex in self.mesh.vertices():
            node = self.mesh.vertex_attribute(vertex, name="node")
            if node is not None:
                xyz = self.graph.node_attributes(node, names="xyz")
                self.mesh.vertex_attributes(vertex, names="xyz", values=xyz)

        angle = radians(self.graph.leaf_angle)

        for loop, vertices in self.mesh.loop_vertices.items():
            points = offset_polyline(self.graph.nodes_attributes(names="xyz", keys=loop), distance=self.node_width)
            for vertex, point in zip(vertices, points):
                self.mesh.vertex_attributes(vertex, names="xyz", values=point)

            if angle != 0:
                u = loop[0]
                uu, vv = vertices[0:2]
                a, b, d = self.mesh.vertices_attributes(names="xyz", keys=[u, uu, vv])
                ab = Line(a, b)
                bd = Line(b, d)
                ab.rotate(angle=-angle, axis=[0, 0, 1], point=a)
                x = intersection_line_line(ab, bd)[0]
                if x:
                    self.mesh.vertex_attributes(uu, names="xyz", values=x)

                u = loop[-1]
                vv, uu = vertices[-2:]
                a, b, d = self.mesh.vertices_attributes(names="xyz", keys=[u, uu, vv])
                ab = Line(a, b)
                bd = Line(b, d)
                ab.rotate(angle=angle, axis=[0, 0, 1], point=a)
                x = intersection_line_line(ab, bd)[0]
                if x:
                    self.mesh.vertex_attributes(uu, names="xyz", values=x)

    def compute_pattern(self) -> None:
        edges = list(flatten(self.mesh.edges_on_boundaries()))
        corners = list(self.mesh.vertices_where(vertex_degree=2))
        self.mesh.edges_attribute(name="crease", value=self.density + 1, keys=edges)
        self._pattern = self.mesh.subdivided(k=self.density, fixed=corners)
