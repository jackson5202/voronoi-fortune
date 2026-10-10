# Voronoi Fortune

Computes the Voronoi diagram of a planar point set using Fortune's sweep line algorithm. Pure JavaScript (ESM), zero dependencies.

## Usage

```js
import { computeVoronoi } from 'voronoi-fortune';

const result = computeVoronoi([
  { x: 0, y: 0 },
  { x: 4, y: 0 },
  { x: 2, y: 3 },
]);

for (const edge of result.edges) {
  console.log(edge.start, edge.end, edge.leftSite, edge.rightSite);
}
```

Pass an explicit clipping box with `{ bbox: { left, right, bottom, top } }`; if omitted, one is derived from the site bounds with padding.

## Why this exists

The standard library offers nothing for computational geometry, and the environments these diagrams are sometimes needed in (locked-down CI, air-gapped containers) forbid installing a native binding. This library is a self-contained implementation small enough to audit and dependency-free enough to run anywhere Node runs.

The trade-off: the beach line is maintained as a linked list rather than a balanced search tree. This keeps the code legible and verifiable, but makes worst-case time O(n²) instead of O(n log n). For the input sizes this library targets — diagrams drawn on screen or embedded in a build step — that is acceptable.

## The awkward edge

Three nearly-collinear sites produce a circumscribed circle whose center is numerically unstable. The implementation uses double precision throughout and tolerates the resulting drift, but if your point set contains points that are almost but not exactly collinear, expect vertex positions to jitter by small amounts. If you need exact arithmetic, this is not the right tool.

Duplicate sites are accepted but produce no edge between them; they collapse into a single cell region.

## Performance

The window keeps a bounded buffer, so `push` is constant time and memory does not
grow with the length of the stream. `peak` and `trough` are linear in the window
size, which is the trade that keeps `push` cheap.

## Limitations

Values are coerced to floats, so very large integers lose precision. If you need
exact integer aggregates over a window, this is the wrong tool.

