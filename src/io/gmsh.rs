//! Gmsh I/O

use crate::traits::{Builder, Entity, Geometry, GmshExport, GmshImport, Grid, Point};
use itertools::izip;
use ndelement::types::ReferenceCellType;
use num::Zero;
use std::collections::HashMap;
use std::fmt::Debug;
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

impl<T: FromStr, B: Builder<T = T, EntityDescriptor = ReferenceCellType>> GmshImport for B
where
    <T as FromStr>::Err: Debug,
{
    fn import_from_gmsh_string(&mut self, s: String) {
        // --- Parse format line ---
        println!("1");
        let format = gmsh_section(&s, "MeshFormat");
        let format_line = format.lines().next().expect("Empty MeshFormat section");
        let mut fmt_parts = format_line.split_whitespace();
        let version = fmt_parts.next().expect("Missing version");
        let ascii_mode = fmt_parts.next().expect("Missing ascii mode");
        println!("2");
        if version != "2.2" && version != "2.0" {
            unimplemented!("Only Gmsh format 2 supported");
        }
        if ascii_mode != "0" {
            unimplemented!("Only ASCII gmsh files supported");
        }
        println!("3");
        // --- Parse nodes ---
        let nodes = gmsh_section(&s, "Nodes");
        let mut lines = nodes.lines();
        let num_nodes: usize = lines
            .next()
            .expect("Nodes section missing num_nodes line")
            .trim()
            .parse()
            .expect("Could not parse num_nodes");
        println!("4");
        for _ in 0..num_nodes {
            let line = lines.next().expect("Unexpected EOF reading nodes");
            let mut parts = line.split_whitespace();
            let tag = parts
                .next()
                .expect("Missing node tag")
                .parse::<usize>()
                .expect("Invalid node tag");
            let coords = parts
                .map(|x| T::from_str(x).expect("Invalid coordinate"))
                .collect::<Vec<_>>();
            if coords.len() != 3 {
                panic!("Expected 3 coordinates per node");
            }
            self.add_point(tag, &coords);
        }
        println!("5");
        // --- Parse elements ---
        let elements = gmsh_section(&s, "Elements");
        let mut lines = elements.lines();
        let num_elements: usize = lines
            .next()
            .expect("Elements section missing num_elements line")
            .trim()
            .parse()
            .expect("Could not parse num_elements");
        println!("6");
        for _ in 0..num_elements {
            let line = lines.next().expect("Unexpected EOF reading elements");
            let parts = line.split_whitespace().collect::<Vec<_>>();

            // Gmsh format 2 element line:
            // element_tag element_type num_tags <tags> node1 node2 ...
            let element_tag: usize = parts[0].parse().unwrap();
            let element_type: usize = parts[1].parse().unwrap();
            let num_tags: usize = parts[2].parse().unwrap();

            let node_start = 3 + num_tags;
            let nodes_line = &parts[node_start..];

            let (cell_type, degree) = interpret_gmsh_cell(element_type);
            let gmsh_perm = get_permutation_to_gmsh(cell_type, degree);

            let mut cell = vec![0; nodes_line.len()];
            for (i, &perm_i) in gmsh_perm.iter().enumerate() {
                cell[perm_i] = nodes_line[i].parse().unwrap();
            }

            self.add_cell_from_nodes_and_type(element_tag, &cell, cell_type, degree);
        }
        println!("7");
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
