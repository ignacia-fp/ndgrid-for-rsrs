//! Geometry where each entity of a given dimension is represented by the same element
use super::{SingleElementGeometry, SingleElementGeometryBorrowed};
use crate::{
    geometry::{Point, PointIter},
    traits::Geometry,
    types::RealScalar,
};
use ndelement::traits::FiniteElement;
use rlst::{RawAccess, Shape};

/// Geometry of an entity
#[derive(Debug)]
pub struct SingleElementEntityGeometry<'a, T: RealScalar, E: FiniteElement> {
    geometry: &'a SingleElementGeometry<T, E>,
    cell_index: usize,
    sub_entity_dimension: usize,
    sub_entity_index: usize,
}

impl<'a, T: RealScalar, E: FiniteElement> SingleElementEntityGeometry<'a, T, E> {
    /// Create new
    pub fn new(
        geometry: &'a SingleElementGeometry<T, E>,
        cell_index: usize,
        sub_entity_dimension: usize,
        sub_entity_index: usize,
    ) -> Self {
        Self {
            geometry,
            cell_index,
            sub_entity_dimension,
            sub_entity_index,
        }
    }
}

impl<T: RealScalar, E: FiniteElement> Geometry for SingleElementEntityGeometry<'_, T, E> {
    type T = T;
    type Point<'a>
        = Point<'a, T>
    where
        Self: 'a;
    type PointIter<'a>
        = PointIter<'a, T>
    where
        Self: 'a;

    fn points(&self) -> PointIter<'_, T> {
        let gdim = self.geometry.dim();
        let mut pts = vec![];
        for index in self
            .geometry
            .element()
            .entity_closure_dofs(self.sub_entity_dimension, self.sub_entity_index)
            .unwrap()
        {
            let i = self.geometry.cells()[[*index, self.cell_index]];
            pts.push((i, &self.geometry.points().data()[i * gdim..(i + 1) * gdim]))
        }

        PointIter::new(pts)
    }

    fn point_count(&self) -> usize {
        self.geometry.cells().shape()[0]
    }
    fn degree(&self) -> usize {
        self.geometry.element().embedded_superdegree()
    }
}
/// Geometry of an entity with borrowed data
#[derive(Debug)]
pub struct SingleElementEntityGeometryBorrowed<'a, T: RealScalar, E: FiniteElement> {
    geometry: &'a SingleElementGeometryBorrowed<'a, T, E>,
    cell_index: usize,
    sub_entity_dimension: usize,
    sub_entity_index: usize,
}

impl<'a, T: RealScalar, E: FiniteElement> SingleElementEntityGeometryBorrowed<'a, T, E> {
    /// Create new
    pub fn new(
        geometry: &'a SingleElementGeometryBorrowed<'a, T, E>,
        cell_index: usize,
        sub_entity_dimension: usize,
        sub_entity_index: usize,
    ) -> Self {
        Self {
            geometry,
            cell_index,
            sub_entity_dimension,
            sub_entity_index,
        }
    }
}

impl<T: RealScalar, E: FiniteElement> Geometry for SingleElementEntityGeometryBorrowed<'_, T, E> {
    type T = T;
    type Point<'a>
        = Point<'a, T>
    where
        Self: 'a;
    type PointIter<'a>
        = PointIter<'a, T>
    where
        Self: 'a;

    fn points(&self) -> PointIter<'_, T> {
        let gdim = self.geometry.dim();
        let mut pts = vec![];
        for index in self
            .geometry
            .element()
            .entity_closure_dofs(self.sub_entity_dimension, self.sub_entity_index)
            .unwrap()
        {
            let i = self.geometry.cells()[[*index, self.cell_index]];
            pts.push((i, &self.geometry.points().data()[i * gdim..(i + 1) * gdim]))
        }

        PointIter::new(pts)
    }

    fn point_count(&self) -> usize {
        self.geometry.cells().shape()[0]
    }
    fn degree(&self) -> usize {
        self.geometry.element().embedded_superdegree()
    }
}
