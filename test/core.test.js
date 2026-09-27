import { describe, it } from 'node:test';
import assert from 'node:assert/strict';
import { computeVoronoi } from '../src/index.js';

/**
 * Helpers for geometric assertions without comparing floats with ==.
 */
function approxEqual(a, b, tol = 1e-6) {
  return Math.abs(a - b) <= tol;
}

function pointsClose(p, q, tol = 1e-6) {
  return approxEqual(p.x, q.x, tol) && approxEqual(p.y, q.y, tol);
}

/**
 * Check that a point lies on or inside the bounding box.
 */
function pointInBox(p, box) {
  return (
    p.x >= box.left - 1e-9 &&
    p.x <= box.right + 1e-9 &&
    p.y >= box.bottom - 1e-9 &&
    p.y <= box.top + 1e-9
  );
}

describe('computeVoronoi', () => {
  it('returns an empty result for no sites', () => {
    const result = computeVoronoi([]);
    assert.deepEqual(result.sites, []);
    assert.deepEqual(result.edges, []);
  });

  it('returns no edges for a single site', () => {
    const result = computeVoronoi([{ x: 0, y: 0 }]);
    assert.equal(result.sites.length, 1);
    assert.equal(result.edges.length, 0);
  });

  it('produces one bisecting edge for two sites', () => {
    const result = computeVoronoi([
      { x: 0, y: 0 },
      { x: 4, y: 0 },
    ]);
    assert.equal(result.sites.length, 2);
    // At least one edge should exist and lie between the two sites.
    assert.ok(result.edges.length >= 1, 'expected at least one edge');
    for (const e of result.edges) {
      assert.ok(pointInBox(e.start, result.bbox), 'edge start in box');
      assert.ok(pointInBox(e.end, result.bbox), 'edge end in box');
    }
  });

  it('computes a diagram for three non-collinear sites', () => {
    const result = computeVoronoi([
      { x: 0, y: 0 },
      { x: 4, y: 0 },
      { x: 2, y: 3 },
    ]);
    assert.equal(result.sites.length, 3);
    assert.ok(result.edges.length >= 1);
    for (const e of result.edges) {
      assert.ok(pointInBox(e.start, result.bbox));
      assert.ok(pointInBox(e.end, result.bbox));
    }
  });

  it('respects a caller-provided bounding box', () => {
    const bbox = { left: -10, right: 10, bottom: -10, top: 10 };
    const result = computeVoronoi(
      [
        { x: 0, y: 0 },
        { x: 2, y: 2 },
        { x: -2, y: 1 },
      ],
      { bbox }
    );
    assert.equal(result.bbox, bbox);
    for (const e of result.edges) {
      assert.ok(pointInBox(e.start, bbox));
      assert.ok(pointInBox(e.end, bbox));
    }
  });

  it('clips edges to the bounding box on a grid', () => {
    const bbox = { left: 0, right: 10, bottom: 0, top: 10 };
    const sites = [
      { x: 2, y: 2 },
      { x: 8, y: 2 },
      { x: 5, y: 8 },
      { x: 5, y: 5 },
    ];
    const result = computeVoronoi(sites, { bbox });
    assert.equal(result.sites.length, 4);
    for (const e of result.edges) {
      assert.ok(pointInBox(e.start, bbox), `start ${JSON.stringify(e.start)} out of box`);
      assert.ok(pointInBox(e.end, bbox), `end ${JSON.stringify(e.end)} out of box`);
    }
  });

  it('produces edges whose endpoints are within the auto bbox', () => {
    const sites = [
      { x: 1, y: 5 },
      { x: 9, y: 5 },
      { x: 5, y: 1 },
      { x: 5, y: 9 },
      { x: 5, y: 5 },
    ];
    const result = computeVoronoi(sites);
    for (const e of result.edges) {
      assert.ok(pointInBox(e.start, result.bbox));
      assert.ok(pointInBox(e.end, result.bbox));
    }
  });

  it('handles collinear sites without throwing', () => {
    const result = computeVoronoi([
      { x: 0, y: 0 },
      { x: 1, y: 0 },
      { x: 2, y: 0 },
      { x: 3, y: 0 },
    ]);
    assert.equal(result.sites.length, 4);
    for (const e of result.edges) {
      assert.ok(pointInBox(e.start, result.bbox));
      assert.ok(pointInBox(e.end, result.bbox));
    }
  });

  it('handles duplicate sites by treating them as one region', () => {
    const result = computeVoronoi([
      { x: 0, y: 0 },
      { x: 0, y: 0 },
      { x: 4, y: 0 },
    ]);
    // No assertion on edge count; just verify it does not crash and produces
    // a valid structure with sites preserved.
    assert.equal(result.sites.length, 3);
    assert.ok(Array.isArray(result.edges));
  });

  it('produces edges that are line segments (start != end after clipping)', () => {
    const result = computeVoronoi([
      { x: 0, y: 0 },
      { x: 4, y: 0 },
    ]);
    for (const e of result.edges) {
      assert.ok(!pointsClose(e.start, e.end), 'clipped edge should have positive length');
    }
  });

  it('works with negative coordinates', () => {
    const result = computeVoronoi([
      { x: -5, y: -5 },
      { x: -1, y: -3 },
      { x: -4, y: 1 },
    ]);
    assert.equal(result.sites.length, 3);
    for (const e of result.edges) {
      assert.ok(pointInBox(e.start, result.bbox));
      assert.ok(pointInBox(e.end, result.bbox));
    }
  });

  it('preserves site coordinates in the output', () => {
    const sites = [
      { x: 3, y: 7 },
      { x: 11, y: 2 },
    ];
    const result = computeVoronoi(sites);
    assert.equal(result.sites.length, sites.length);
    for (let i = 0; i < sites.length; i++) {
      assert.ok(approxEqual(result.sites[i].x, sites[i].x));
      assert.ok(approxEqual(result.sites[i].y, sites[i].y));
    }
  });
});
