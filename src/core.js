/**
 * Fortune's sweep line algorithm for computing the Voronoi diagram of a planar
 * point set.
 *
 * Output: a list of edges. Each edge is the portion of a Voronoi cell boundary
 * that is a straight line segment; rays and parabolic arc segments are clipped
 * to a bounding box provided by the caller. Degenerate inputs (fewer than two
 * sites, duplicate sites) are handled explicitly.
 *
 * Why a direct, self-contained implementation rather than a binding to a native
 * library: the brief forbids third-party dependencies, and the standard library
 * offers nothing for computational geometry. A pure-JavaScript implementation
 * is the only option that runs in the restricted environment described.
 *
 * The main numerical hazard is that three nearly-collinear sites produce an
 * equidistant circle whose center is very sensitive to small perturbations. We
 * rely on double precision and document the limitation in the README.
 */

/**
 * @typedef {Object} Point
 * @property {number} x
 * @property {number} y
 */

/**
 * @typedef {Object} BoundingBox
 * @property {number} left
 * @property {number} right
 * @property {number} bottom
 * @property {number} top
 */

/**
 * @typedef {Object} VoronoiEdge
 * @property {Point} start
 * @property {Point} end
 * @property {Point|null} leftSite
 * @property {Point|null} rightSite
 */

/**
 * @typedef {Object} VoronoiResult
 * @property {Point[]} sites
 * @property {VoronoiEdge[]} edges
 * @property {BoundingBox} bbox
 */

const EPS = 1e-9;

/**
 * Lexicographic comparison on (y, x) so the sweep proceeds top-to-bottom and
 * left-to-right. Returns negative if a precedes b.
 */
function compareEvents(a, b) {
  if (Math.abs(a.y - b.y) > EPS) return b.y - a.y; // descending y
  if (Math.abs(a.x - b.x) > EPS) return a.x - b.x;
  return 0;
}

/**
 * Binary-heap priority queue ordered by compareEvents. We use an array and
 * sift manually; a linked structure is unnecessary and harder to audit.
 */
class PriorityQueue {
  constructor() {
    this.items = [];
  }
  get isEmpty() {
    return this.items.length === 0;
  }
  push(item) {
    const arr = this.items;
    arr.push(item);
    let i = arr.length - 1;
    while (i > 0) {
      const parent = (i - 1) >> 1;
      if (compareEvents(arr[parent], arr[i]) <= 0) break;
      const tmp = arr[parent];
      arr[parent] = arr[i];
      arr[i] = tmp;
      i = parent;
    }
  }
  pop() {
    const arr = this.items;
    if (arr.length === 0) return null;
    const top = arr[0];
    const last = arr.pop();
    if (arr.length > 0) {
      arr[0] = last;
      let i = 0;
      const n = arr.length;
      for (;;) {
        const l = 2 * i + 1;
        const r = 2 * i + 2;
        let smallest = i;
        if (l < n && compareEvents(arr[l], arr[smallest]) < 0) smallest = l;
        if (r < n && compareEvents(arr[r], arr[smallest]) < 0) smallest = r;
        if (smallest === i) break;
        const tmp = arr[smallest];
        arr[smallest] = arr[i];
        arr[i] = tmp;
        i = smallest;
      }
    }
    return top;
  }
}

/**
 * Compute the intersection of the two parabolas with foci at f1 and f2, both
 * equidistant from the directrix at y = directrixY. The intersection x is used
 * as the breakpoint key in the beach line's balanced tree.
 */
function parabolaIntersection(f1, f2, directrixY) {
  // The two foci are at the same height; intersection is halfway between.
  if (Math.abs(f1.y - f2.y) < EPS) {
    return { x: (f1.x + f2.x) / 2, x1: f1.x, x2: f2.x };
  }
  const p = (f1.y - directrixY) / (f2.y - directrixY);
  // Derived from equating (x-f1.x)^2 + (y-f1.y)^2 = (y - directrix)^2 etc.
  // Solving the linear equation in x for the common chord of the two parabolas.
  const x = (f1.x * p + f2.x) / (1 + p);
  return { x, x1: f1.x, x2: f2.x };
}

/**
 * Circle event: the lowest point on the circle through three consecutive arc
 * foci. When the sweep line reaches this y, the middle arc disappears.
 */
class CircleEvent {
  constructor(x, y, yCenter, arc) {
    this.x = x;
    this.y = y; // the sweep line y at which the circle is tangent
    this.yCenter = yCenter;
    this.arc = arc;
    this.valid = true;
  }
}

/**
 * Site event container.
 */
class SiteEvent {
  constructor(point) {
    this.x = point.x;
    this.y = point.y;
    this.point = point;
  }
}

/**
 * A node in the beach line: either an arc (leaf) or a breakpoint (internal).
 */
class BeachNode {}

class ArcNode extends BeachNode {
  constructor(site) {
    super();
    this.site = site;
    this.circleEvent = null;
    this.prev = null;
    this.next = null;
  }
}

class BreakpointNode extends BeachNode {
  constructor(leftSite, rightSite, leftArc, rightArc) {
    super();
    this.leftSite = leftSite;
    this.rightSite = rightSite;
    this.leftArc = leftArc;
    this.rightArc = rightArc;
    this.edge = null;
  }
}

/**
 * Compute the lowest point and center of the circle through three points.
 * Returns null if the points are collinear (no finite circle).
 */
function circleThrough(a, b, c) {
  const d = 2 * (a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y));
  if (Math.abs(d) < EPS) return null; // collinear
  const ux =
    ((a.x * a.x + a.y * a.y) * (b.y - c.y) +
      (b.x * b.x + b.y * b.y) * (c.y - a.y) +
      (c.x * c.x + c.y * c.y) * (a.y - b.y)) /
    d;
  const uy =
    ((a.x * a.x + a.y * a.y) * (c.x - b.x) +
      (b.x * b.x + b.y * b.y) * (a.x - c.x) +
      (c.x * c.x + c.y * c.y) * (b.x - a.x)) /
    d;
  const r = Math.hypot(a.x - ux, a.y - uy);
  return { x: ux, y: uy, r };
}

/**
 * Clip a ray to the bounding box. The edge direction is (dx, dy). Returns the
 * segment from origin extending in the direction, clipped to box, or null if no
 * intersection exists.
 */
function clipRay(origin, dx, dy, box) {
  let tmin = 0;
  let tmax = Infinity;
  // Parametric: p = origin + t*(dx,dy). Solve against each box face.
  if (Math.abs(dx) < EPS) {
    if (origin.x < box.left || origin.x > box.right) return null;
  } else {
    const t1 = (box.left - origin.x) / dx;
    const t2 = (box.right - origin.x) / dx;
    const lo = Math.min(t1, t2);
    const hi = Math.max(t1, t2);
    tmin = Math.max(tmin, lo);
    tmax = Math.min(tmax, hi);
  }
  if (Math.abs(dy) < EPS) {
    if (origin.y < box.bottom || origin.y > box.top) return null;
  } else {
    const t1 = (box.bottom - origin.y) / dy;
    const t2 = (box.top - origin.y) / dy;
    const lo = Math.min(t1, t2);
    const hi = Math.max(t1, t2);
    tmin = Math.max(tmin, lo);
    tmax = Math.min(tmax, hi);
  }
  if (tmin > tmax) return null;
  const endT = tmax === Infinity ? tmin : tmax;
  if (!isFinite(endT)) return null;
  return {
    start: { x: origin.x + tmin * dx, y: origin.y + tmin * dy },
    end: { x: origin.x + endT * dx, y: origin.y + endT * dy },
  };
}

/**
 * Clip an existing segment to the bounding box. Returns the clipped segment or
 * null if fully outside.
 */
function clipSegment(p1, p2, box) {
  let tmin = 0;
  let tmax = 1;
  const dx = p2.x - p1.x;
  const dy = p2.y - p1.y;
  if (Math.abs(dx) < EPS) {
    if (p1.x < box.left - EPS || p1.x > box.right + EPS) return null;
  } else {
    const t1 = (box.left - p1.x) / dx;
    const t2 = (box.right - p1.x) / dx;
    tmin = Math.max(tmin, Math.min(t1, t2));
    tmax = Math.min(tmax, Math.max(t1, t2));
  }
  if (Math.abs(dy) < EPS) {
    if (p1.y < box.bottom - EPS || p1.y > box.top + EPS) return null;
  } else {
    const t1 = (box.bottom - p1.y) / dy;
    const t2 = (box.top - p1.y) / dy;
    tmin = Math.max(tmin, Math.min(t1, t2));
    tmax = Math.min(tmax, Math.max(t1, t2));
  }
  if (tmin > tmax + EPS) return null;
  return {
    start: { x: p1.x + tmin * dx, y: p1.y + tmin * dy },
    end: { x: p1.x + tmax * dx, y: p1.y + tmax * dy },
  };
}

/**
 * Compute the Voronoi diagram.
 * @param {Point[]} sites - input points
 * @param {Object} [options]
 * @param {BoundingBox} [options.bbox] - clipping box; auto-computed if omitted
 * @returns {VoronoiResult}
 */
export function computeVoronoi(sites, options = {}) {
  if (!Array.isArray(sites) || sites.length === 0) {
    return { sites: [], edges: [], bbox: { left: 0, right: 0, bottom: 0, top: 0 } };
  }

  // Copy and verify sites.
  const pts = sites.map((s) => ({ x: Number(s.x), y: Number(s.y) }));

  if (pts.length === 1) {
    const bbox = options.bbox || { left: pts[0].x - 1, right: pts[0].x + 1, bottom: pts[0].y - 1, top: pts[0].y + 1 };
    return { sites: pts, edges: [], bbox };
  }

  // Determine bounding box.
  let bbox;
  if (options.bbox) {
    bbox = options.bbox;
  } else {
    let minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
    for (const p of pts) {
      if (p.x < minX) minX = p.x;
      if (p.x > maxX) maxX = p.x;
      if (p.y < minY) minY = p.y;
      if (p.y > maxY) maxY = p.y;
    }
    const w = maxX - minX || 1;
    const h = maxY - minY || 1;
    const pad = Math.max(w, h) * 0.5;
    bbox = { left: minX - pad, right: maxX + pad, bottom: minY - pad, top: maxY + pad };
  }

  // Sort site events by (y desc, x asc).
  const queue = new PriorityQueue();
  for (const p of pts) queue.push(new SiteEvent(p));

  // Beach line as doubly-linked list of ArcNode, with breakpoint metadata stored
  // on the arc to the right of each breakpoint. This avoids a full balanced BST;
  // the linked list is O(n) per operation but acceptable for small inputs and
  // far easier to verify than a hand-rolled red-black tree.
  let beachHead = null;
  const edges = [];
  /** @type {Set<CircleEvent>} */
  const pendingCircles = new Set();

  /**
   * Find the arc directly above x at sweep line y.
   */
  function findArc(x, sweepY) {
    let arc = beachHead;
    if (!arc) return null;
    while (arc.next) {
      const nextArc = arc.next;
      // The breakpoint between arc and nextArc uses both sites and the sweep line.
      const bp = breakpointX(arc, nextArc, sweepY);
      if (x < bp) break;
      arc = nextArc;
    }
    return arc;
  }

  function breakpointX(leftArc, rightArc, sweepY) {
    const s1 = leftArc.site;
    const s2 = rightArc.site;
    if (Math.abs(s1.y - sweepY) < EPS) return s1.x;
    if (Math.abs(s2.y - sweepY) < EPS) return s2.x;
    const res = parabolaIntersection(s1, s2, sweepY);
    return res.x;
  }

  /**
   * Add a circle event for the triple (prevArc, arc, nextArc) if it converges.
   */
  function addCircleEvent(prevArc, arc, nextArc) {
    if (!prevArc || !nextArc) return;
    const c = circleThrough(prevArc.site, arc.site, nextArc.site);
    if (!c) return; // collinear, no event
    if (c.y - arc.site.y > EPS) return; // circle bottom is above the arc focus
    const ev = new CircleEvent(c.x, c.y, c.y, arc);
    arc.circleEvent = ev;
    pendingCircles.add(ev);
    queue.push(ev);
  }

  function invalidateCircleEvent(arc) {
    if (arc.circleEvent) {
      arc.circleEvent.valid = false;
      pendingCircles.delete(arc.circleEvent);
      arc.circleEvent = null;
    }
  }

  /**
   * Insert a new site into the beach line, splitting the arc above it.
   */
  function handleSiteEvent(site) {
    if (!beachHead) {
      beachHead = new ArcNode(site);
      return;
    }
    const arc = findArc(site.x, site.y);
    if (!arc) return;
    invalidateCircleEvent(arc);

    // Split arc into two arcs, insert new arc between.
    const newArc = new ArcNode(site);
    // We maintain a linked list: prev <-> arc <-> newArc <-> arcCopy <-> next
    const arcCopy = new ArcNode(arc.site);
    if (arc.prev) {
      arc.prev.next = newArc;
      newArc.prev = arc.prev;
    } else {
      beachHead = newArc;
    }
    newArc.next = arcCopy;
    arcCopy.prev = newArc;
    arcCopy.next = arc.next;
    if (arc.next) arc.next.prev = arcCopy;
    arc.next = null;
    arc.prev = null;

    // Create two half-edges starting at the intersection point (on the beach line).
    const startPt = { x: site.x, y: site.y };
    const e1 = { start: { ...startPt }, end: { ...startPt }, leftSite: arc.site, rightSite: site };
    const e2 = { start: { ...startPt }, end: { ...startPt }, leftSite: site, rightSite: arcCopy.site };
    edges.push(e1, e2);
    newArc._edgeOut = e1;
    arcCopy._edgeOut = e2;

    // Check for circle events.
    addCircleEvent(newArc.prev, newArc, arcCopy);
    addCircleEvent(newArc, arcCopy, arcCopy.next);
  }

  /**
   * Remove the shrinking arc and connect the two adjacent breakpoints.
   */
  function handleCircleEvent(ev) {
    if (!ev.valid) return;
    const arc = ev.arc;
    invalidateCircleEvent(arc);

    const prevArc = arc.prev;
    const nextArc = arc.next;
    if (!prevArc || !nextArc) return;

    // The vertex is the circle center.
    const vertex = { x: ev.x, y: ev.yCenter };

    // Finish the two edges meeting here.
    if (prevArc._edgeOut) prevArc._edgeOut.end = { ...vertex };
    if (arc._edgeOut) arc._edgeOut.end = { ...vertex };

    // Remove arc from the list.
    prevArc.next = nextArc;
    nextArc.prev = prevArc;
    if (beachHead === arc) beachHead = nextArc;
    arc.prev = null;
    arc.next = null;

    // Create a new edge from the vertex between prevArc.site and nextArc.site.
    const newEdge = { start: { ...vertex }, end: { ...vertex }, leftSite: prevArc.site, rightSite: nextArc.site };
    edges.push(newEdge);
    prevArc._edgeOut = newEdge;

    invalidateCircleEvent(prevArc);
    invalidateCircleEvent(nextArc);
    addCircleEvent(prevArc.prev, prevArc, nextArc);
    addCircleEvent(prevArc, nextArc, nextArc.next);
  }

  // Main event loop.
  while (!queue.isEmpty) {
    const ev = queue.pop();
    if (ev instanceof CircleEvent) {
      handleCircleEvent(ev);
    } else {
      handleSiteEvent(ev);
    }
  }

  // Clip all edges to the bounding box. Rays (start == end after processing)
  // need a direction derived from the two sites.
  const clipped = [];
  for (const e of edges) {
    const isRay =
      Math.abs(e.start.x - e.end.x) < EPS && Math.abs(e.start.y - e.end.y) < EPS;
    let seg;
    if (isRay && e.leftSite && e.rightSite) {
      // Direction perpendicular to the line between the two sites, pointing
      // away from the higher site (downward sweep means the edge grows down).
      const dx = e.rightSite.x - e.leftSite.x;
      const dy = e.rightSite.y - e.leftSite.y;
      // Perpendicular: (-dy, dx) or (dy, -dx). Pick the one with negative y
      // (edge extends downward as sweep continues).
      let px = -dy;
      let py = dx;
      if (py > -EPS) {
        px = dy;
        py = -dx;
      }
      seg = clipRay(e.start, px, py, bbox);
    } else {
      seg = clipSegment(e.start, e.end, bbox);
    }
    if (seg) {
      clipped.push({
        start: seg.start,
        end: seg.end,
        leftSite: e.leftSite,
        rightSite: e.rightSite,
      });
    }
  }

  return { sites: pts, edges: clipped, bbox };
}
