//! Gmsh I/O

use crate::traits::{Builder, Entity, Geometry, GmshExport, GmshImport, Grid, Point};
use itertools::izip;
use ndelement::types::ReferenceCellType;
use num::Zero;
use std::collections::HashMap;
use std::str::FromStr;

fn get_permutation_to_gmsh(cell_type: ReferenceCellType, degree: usize) -> Vec<usize> {
    match cell_type {
        ReferenceCellType::Triangle => match degree {
            1 => vec![0, 1, 2],
            2 => vec![0, 1, 2, 5, 3, 4],
            3 => vec![0, 1, 2, 7, 8, 3, 4, 6, 5, 9],
            4 => vec![0, 1, 2, 9, 10, 11, 3, 4, 5, 8, 7, 6, 12, 13, 14],
            5 => vec![
                0, 1, 2, 11, 12, 13, 14, 3, 4, 5, 6, 10, 9, 8, 7, 15, 16, 17, 18, 19, 20,
            ],
            _ => {
                panic!("Unsupported degree");
            }
        },
        ReferenceCellType::Quadrilateral => match degree {
            1 => vec![0, 1, 3, 2],
            2 => vec![0, 1, 3, 2, 4, 6, 7, 5, 8],
            _ => {
                panic!("Unsupported degree");
            }
        },
        ReferenceCellType::Tetrahedron => match degree {
            1 => vec![0, 1, 2, 3],
            _ => {
                panic!("Unsupported degree");
            }
        },
        ReferenceCellType::Hexahedron => match degree {
            1 => vec![0, 1, 3, 2, 4, 5, 7, 6],
            _ => {
                panic!("Unsupported degree");
            }
        },
        _ => {
            panic!("Unsupported cell type.");
        }
    }
}

fn get_gmsh_cell(cell_type: ReferenceCellType, degree: usize) -> usize {
    match cell_type {
        ReferenceCellType::Triangle => match degree {
            1 => 2,
            2 => 9,
            3 => 21,
            4 => 23,
            5 => 25,
            _ => {
                panic!("Unsupported degree");
            }
        },
        ReferenceCellType::Quadrilateral => match degree {
            1 => 3,
            2 => 10,
            _ => {
                panic!("Unsupported degree");
            }
        },
        ReferenceCellType::Tetrahedron => match degree {
            1 => 4,
            _ => {
                panic!("Unsupported degree");
            }
        },
        ReferenceCellType::Hexahedron => match degree {
            1 => 5,
            _ => {
                panic!("Unsupported degree");
            }
        },
        _ => {
            panic!("Unsupported cell type.");
        }
    }
}

fn interpret_gmsh_cell(gmsh_cell: usize) -> (ReferenceCellType, usize) {
    match gmsh_cell {
        2 => (ReferenceCellType::Triangle, 1),
        9 => (ReferenceCellType::Triangle, 2),
        21 => (ReferenceCellType::Triangle, 3),
        23 => (ReferenceCellType::Triangle, 4),
        25 => (ReferenceCellType::Triangle, 5),
        3 => (ReferenceCellType::Quadrilateral, 1),
        10 => (ReferenceCellType::Quadrilateral, 2),
        4 => (ReferenceCellType::Tetrahedron, 1),
        5 => (ReferenceCellType::Hexahedron, 1),
        _ => {
            panic!("Unsupported cell type.");
        }
    }
}

impl<G: Grid<EntityDescriptor = ReferenceCellType>> GmshExport for G {
    fn to_gmsh_string(&self) -> String {
        let tdim = self.topology_dim();
        let gdim = self.geometry_dim();

        let mut points = HashMap::new();
        for cell in self.cell_iter() {
            for point in cell.geometry().points() {
                let mut p = vec![G::T::zero(); gdim];
                point.coords(&mut p);
                points.insert(point.index(), p);
            }
        }
        let mut points = points.iter().collect::<Vec<_>>();
        points.sort_by(|i, j| i.0.cmp(j.0));
        let node_count = points.len();

        let mut gmsh_s = String::from("");
        gmsh_s.push_str("$MeshFormat\n");
        gmsh_s.push_str("4.1 0 8\n");
        gmsh_s.push_str("$EndMeshFormat\n");
        gmsh_s.push_str("$Nodes\n");
        gmsh_s.push_str(&format!("1 {node_count} 1 {node_count}\n"));
        gmsh_s.push_str(&format!("{tdim} 1 0 {node_count}\n"));
        for (i, _) in &points {
            gmsh_s.push_str(&format!("{}\n", *i + 1));
        }
        for (_, coords) in &points {
            for (n, c) in coords.iter().enumerate() {
                if n != 0 {
                    gmsh_s.push(' ');
                }
                gmsh_s.push_str(&format!("{c}"));
            }
            gmsh_s.push('\n');
        }
        gmsh_s.push_str("$EndNodes\n");
        gmsh_s.push_str("$Elements\n");

        let cell_count = self
            .entity_types(tdim)
            .iter()
            .map(|t| self.entity_count(*t))
            .sum::<usize>();

        let mut elements = vec![];
        let mut cells_by_element = vec![];
        for (i, cell) in self.entity_iter(tdim).enumerate() {
            let element = (cell.entity_type(), cell.geometry().degree());
            if !elements.contains(&element) {
                elements.push(element);
                cells_by_element.push(vec![]);
            }
            cells_by_element[elements.iter().position(|i| *i == element).unwrap()].push(i);
        }

        gmsh_s.push_str(&format!("{} {cell_count} 1 {cell_count}\n", elements.len()));

        for ((cell_type, degree), cells) in izip!(elements, cells_by_element) {
            let gmsh_perm = get_permutation_to_gmsh(cell_type, degree);

            gmsh_s.push_str(&format!(
                "{tdim} 1 {} {}\n",
                get_gmsh_cell(cell_type, degree),
                cells.len()
            ));

            for (i, index) in cells.iter().enumerate() {
                gmsh_s.push_str(&format!("{}", i + 1));
                let entity = self.entity(tdim, *index).unwrap();
                let point_indices = entity
                    .geometry()
                    .points()
                    .map(|i| i.index())
                    .collect::<Vec<_>>();
                for j in &gmsh_perm {
                    gmsh_s.push_str(&format!(" {}", point_indices[*j] + 1));
                }
                gmsh_s.push('\n');
            }
        }

        gmsh_s.push_str("$EndElements\n");

        gmsh_s
    }
}

/// Get a section from a gmsh string
fn gmsh_section(s: &str, section: &str) -> String {
    let a = s.split(&format!("${section}\n")).collect::<Vec<_>>();
    if a.len() <= 1 {
        panic!("Section not found: {section}");
    }
    String::from(a[1].split(&format!("\n$End{section}")).collect::<Vec<_>>()[0])
}

impl<
        T: FromStr + Zero + std::marker::Copy,
        B: Builder<T = T, EntityDescriptor = ReferenceCellType>,
    > GmshImport for B
where
    <T as FromStr>::Err: std::fmt::Debug,
{
    fn import_from_gmsh_string(&mut self, s: String) {
        // --- Parse MeshFormat section ---
        let format = gmsh_section(&s, "MeshFormat");
        let mut format_parts = format.split_whitespace();
        let version = format_parts.next().unwrap();
        let ascii_mode = format_parts.next().unwrap();

        if version != "4.1" {
            unimplemented!("Unsupported gmsh file version");
        }
        if ascii_mode != "0" {
            unimplemented!("Non-ASCII gmsh files currently not supported");
        }

        // --- Parse Nodes section ---
        let nodes = gmsh_section(&s, "Nodes");
        let nodes: Vec<&str> = nodes.lines().collect();

        let mut header_parts = nodes[0].split_whitespace();
        let num_entity_blocks = header_parts.next().unwrap().parse::<usize>().unwrap();
        // _num_nodes, _min_node_tag, _max_node_tag are unused

        let mut line_n = 1;
        let mut coords = [T::zero(); 3];
        for _ in 0..num_entity_blocks {
            let mut block_header = nodes[line_n].split_whitespace();
            let _entity_dim = block_header.next().unwrap().parse::<usize>().unwrap();
            let _entity_tag = block_header.next().unwrap().parse::<usize>().unwrap();
            let parametric = block_header.next().unwrap().parse::<usize>().unwrap();
            let num_nodes_in_block = block_header.next().unwrap().parse::<usize>().unwrap();

            if parametric == 1 {
                unimplemented!("Parametric nodes currently not supported");
            }

            line_n += 1;
            let tags = &nodes[line_n..line_n + num_nodes_in_block];
            let coords_lines = &nodes[line_n + num_nodes_in_block..line_n + 2 * num_nodes_in_block];

            for (t_line, c_line) in izip!(tags, coords_lines) {
                let tag = t_line.trim().parse::<usize>().unwrap();

                for (i, val) in c_line.trim().split_whitespace().enumerate() {
                    coords[i] = T::from_str(val).expect("Could not parse coordinate");
                }

                self.add_point(tag, &coords);
            }

            line_n += num_nodes_in_block * 2;
        }

        // --- Parse Elements section ---
        let elements = gmsh_section(&s, "Elements");
        let elements: Vec<&str> = elements.lines().collect();

        let mut header_parts = elements[0].split_whitespace();
        let num_entity_blocks = header_parts.next().unwrap().parse::<usize>().unwrap();
        // _num_elements, _min_element_tag, _max_element_tag are unused

        line_n = 1;
        let mut cell = vec![0usize; 32]; // Arbitrary max size, will be trimmed

        for _ in 0..num_entity_blocks {
            let mut block_header = elements[line_n].split_whitespace();
            let _entity_dim = block_header.next().unwrap().parse::<usize>().unwrap();
            let _entity_tag = block_header.next().unwrap().parse::<usize>().unwrap();
            let element_type = block_header.next().unwrap().parse::<usize>().unwrap();
            let num_elements_in_block = block_header.next().unwrap().parse::<usize>().unwrap();

            let (cell_type, degree) = interpret_gmsh_cell(element_type);
            let gmsh_perm = get_permutation_to_gmsh(cell_type, degree);
            let num_nodes_per_cell = gmsh_perm.len();

            line_n += 1;
            for line in &elements[line_n..line_n + num_elements_in_block] {
                let mut iter = line.split_whitespace().map(|i| i.parse::<usize>().unwrap());
                let tag = iter.next().unwrap();

                let nodes: Vec<usize> = iter.collect(); // usually small, e.g. 3 or 4
                debug_assert_eq!(nodes.len(), num_nodes_per_cell);

                for (i, &j) in gmsh_perm.iter().enumerate() {
                    cell[i] = nodes[j];
                }

                self.add_cell_from_nodes_and_type(
                    tag,
                    &cell[..num_nodes_per_cell],
                    cell_type,
                    degree,
                );
            }

            line_n += num_elements_in_block;
        }
    }

    fn import_from_gmsh(&mut self, filename: &str) {
        let content = std::fs::read_to_string(filename).expect("Unable to read file");
        self.import_from_gmsh_string(content);
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::{
        shapes::regular_sphere,
        traits::{Builder, Topology},
        SingleElementGridBuilder,
    };
    use approx::*;

    #[test]
    fn test_regular_sphere_gmsh_io() {
        let g = regular_sphere::<f64>(2);
        g.export_as_gmsh("_test_io_sphere.msh");
    }

    #[test]
    fn test_export_quads() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Quadrilateral, 1));
        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[1.0, 0.0, 0.0]);
        b.add_point(2, &[0.0, 1.0, 0.0]);
        b.add_point(3, &[1.0, 1.0, 0.0]);
        b.add_point(4, &[0.0, 0.0, 1.0]);
        b.add_point(5, &[1.0, 0.0, 1.0]);
        b.add_point(6, &[0.0, 1.0, 1.0]);
        b.add_point(7, &[1.0, 1.0, 1.0]);
        b.add_cell(0, &[0, 2, 1, 3]);
        b.add_cell(1, &[0, 1, 4, 5]);
        b.add_cell(2, &[0, 4, 2, 6]);
        b.add_cell(3, &[1, 3, 5, 7]);
        b.add_cell(4, &[2, 6, 3, 7]);
        b.add_cell(5, &[4, 5, 6, 7]);
        let g = b.create_grid();
        g.export_as_gmsh("_test_io_cube.msh");
    }

    #[test]
    fn test_export_tetrahedra() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Tetrahedron, 1));
        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[1.0, 0.0, 0.0]);
        b.add_point(2, &[0.0, 1.0, 0.0]);
        b.add_point(3, &[1.0, 1.0, 0.0]);
        b.add_point(4, &[0.0, 0.0, 1.0]);
        b.add_point(5, &[1.0, 0.0, 1.0]);
        b.add_point(6, &[0.0, 1.0, 1.0]);
        b.add_point(7, &[1.0, 1.0, 1.0]);
        b.add_cell(0, &[0, 1, 5, 7]);
        b.add_cell(1, &[0, 2, 6, 7]);
        b.add_cell(2, &[0, 4, 5, 7]);
        b.add_cell(3, &[0, 1, 3, 7]);
        b.add_cell(4, &[0, 2, 3, 7]);
        b.add_cell(5, &[0, 4, 6, 7]);
        let g = b.create_grid();
        g.export_as_gmsh("_test_io_tetrahedra.msh");
    }

    #[test]
    fn test_hexahedra() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Hexahedron, 1));
        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[1.0, 0.0, 0.0]);
        b.add_point(2, &[0.0, 1.0, 0.0]);
        b.add_point(3, &[1.0, 1.0, 0.0]);
        b.add_point(4, &[0.0, 0.0, 1.0]);
        b.add_point(5, &[1.0, 0.0, 1.0]);
        b.add_point(6, &[0.0, 1.0, 1.0]);
        b.add_point(7, &[1.0, 1.0, 1.0]);
        b.add_point(8, &[0.0, 0.0, 2.0]);
        b.add_point(9, &[1.0, 0.0, 2.0]);
        b.add_point(10, &[0.0, 1.0, 2.0]);
        b.add_point(11, &[1.0, 1.0, 1.5]);
        b.add_cell(1, &[0, 1, 2, 3, 4, 5, 6, 7]);
        b.add_cell(2, &[4, 5, 6, 7, 8, 9, 10, 11]);
        let g = b.create_grid();
        g.export_as_gmsh("_test_io_hexahedra.msh");

        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Hexahedron, 1));
        b.import_from_gmsh("_test_io_hexahedra.msh");
        let g2 = b.create_grid();

        let mut p1 = [0.0; 3];
        let mut p2 = [0.0; 3];
        for (v1, v2) in izip!(g.entity_iter(0), g2.entity_iter(0)) {
            for (pt1, pt2) in izip!(v1.geometry().points(), v2.geometry().points()) {
                pt1.coords(&mut p1);
                pt2.coords(&mut p2);
                for (c1, c2) in izip!(&p1, &p2) {
                    assert_relative_eq!(c1, c2, epsilon = 1e-10);
                }
            }
        }
        for (h1, h2) in izip!(g.entity_iter(3), g2.entity_iter(3)) {
            for (v1, v2) in izip!(
                h1.topology().sub_entity_iter(0),
                h2.topology().sub_entity_iter(0)
            ) {
                assert_eq!(v1, v2);
            }
        }
    }

    #[test]
    fn test_high_order_triangles() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 5));
        b.add_point(0, &[0.0, 0.0, 0.0]);
        b.add_point(1, &[5.0, 0.0, 0.2]);
        b.add_point(2, &[0.0, 5.0, 0.4]);
        b.add_point(3, &[4.0, 1.0, -0.1]);
        b.add_point(4, &[3.0, 2.0, 0.2]);
        b.add_point(5, &[2.0, 3.0, -0.3]);
        b.add_point(6, &[1.0, 4.0, -0.5]);
        b.add_point(7, &[0.0, 1.0, 0.6]);
        b.add_point(8, &[0.0, 2.0, 0.2]);
        b.add_point(9, &[0.0, 3.0, 0.1]);
        b.add_point(10, &[0.0, 4.0, -0.2]);
        b.add_point(11, &[1.0, 0.0, -0.3]);
        b.add_point(12, &[2.0, 0.0, -0.4]);
        b.add_point(13, &[3.0, 0.0, -0.5]);
        b.add_point(14, &[4.0, 0.0, -0.2]);
        b.add_point(15, &[1.0, 1.0, 0.1]);
        b.add_point(16, &[2.0, 1.0, 0.1]);
        b.add_point(17, &[3.0, 1.0, 0.1]);
        b.add_point(18, &[2.0, 1.0, 0.2]);
        b.add_point(19, &[2.0, 2.0, 0.1]);
        b.add_point(20, &[3.0, 1.0, 0.1]);
        b.add_cell(
            1,
            &[
                0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20,
            ],
        );
        let g = b.create_grid();
        g.export_as_gmsh("_test_io_high_order_triangle.msh");

        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 5));
        b.import_from_gmsh("_test_io_high_order_triangle.msh");
        let g2 = b.create_grid();

        let mut p1 = [0.0; 3];
        let mut p2 = [0.0; 3];
        for (v1, v2) in izip!(g.entity_iter(0), g2.entity_iter(0)) {
            v1.geometry().points().next().unwrap().coords(&mut p1);
            v2.geometry().points().next().unwrap().coords(&mut p2);
        }
        for (v1, v2) in izip!(g.entity_iter(0), g2.entity_iter(0)) {
            v1.geometry().points().next().unwrap().coords(&mut p1);
            v2.geometry().points().next().unwrap().coords(&mut p2);
            for (c1, c2) in izip!(&p1, &p2) {
                assert_relative_eq!(c1, c2, epsilon = 1e-10);
            }
        }
        for (h1, h2) in izip!(g.entity_iter(2), g2.entity_iter(2)) {
            for (v1, v2) in izip!(
                h1.topology().sub_entity_iter(0),
                h2.topology().sub_entity_iter(0)
            ) {
                assert_eq!(v1, v2);
            }
        }
    }

    #[test]
    fn test_import_triangle() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 1));
        b.import_from_gmsh("meshes/sphere_triangle.msh");
        let _g = b.create_grid();
    }

    #[test]
    fn test_import_quadrilateral() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Quadrilateral, 1));
        b.import_from_gmsh("meshes/cube_quadrilateral.msh");
        let _g = b.create_grid();
    }

    #[test]
    fn test_import_tetrahedron() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Tetrahedron, 1));
        b.import_from_gmsh("meshes/cube_tetrahedron.msh");
        let _g = b.create_grid();
    }

    #[test]
    #[should_panic]
    fn test_import_wrong_cell() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Quadrilateral, 1));
        b.import_from_gmsh("meshes/cube_tetrahedron.msh");
        let _g = b.create_grid();
    }

    #[test]
    #[should_panic]
    fn test_import_wrong_degree() {
        let mut b = SingleElementGridBuilder::<f64>::new(3, (ReferenceCellType::Triangle, 2));
        b.import_from_gmsh("meshes/sphere_triangle.msh");
        let _g = b.create_grid();
    }
}
