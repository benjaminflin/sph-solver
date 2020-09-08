import Vector2 from "./Vector2";

const canvas = <HTMLCanvasElement>document.getElementById("canvas");
const ctx = canvas.getContext("2d");

canvas.width = window.innerWidth - 10;
canvas.height = window.innerHeight - 10;

const REST_DENS = 1000; // rest density
const GAS_CONST = 2000; // const for equation of state
const H = 16; // kernel radius
const HSQ = H * H; // radius^2 for optimization
const MASS = 65; // assume all particles have the same mass
const VISC = 250; // viscosity constant
const DT = 0.0008; // integration timestep
const BOUND_DAMPING = -0.5;
const POLY6 = 315 / (65 * Math.PI * Math.pow(H, 9));
// const SPIKY_GRAD = 0.5;
const SPIKY_GRAD = -45.0 / (Math.PI * Math.pow(H, 6));
const VISC_LAP = 45 / (Math.PI * Math.pow(H, 6));
const SIGMA = 100000000;
// particles size
const P_SIZE = 5;
const G = new Vector2(0, 12000 * 9.8);
// number of particles
const NUM_PARTICLES = 1000;
const EPSILON = 0.001;
// particles
let particles: Particle[] = [];

let grid: Particle[][] = [];
let num_rows = 0;
let num_cols = 0;

class Particle {
  position: Vector2;
  velocity: Vector2;
  force: Vector2;
  density: number;
  pressure: number;
  constructor(position = new Vector2()) {
    this.position = position;
    this.velocity = new Vector2();
    this.force = new Vector2();
    this.density = 0;
    this.pressure = 0;
  }
}

const initGrid = () => {
  num_rows = Math.ceil(canvas.width / (2 * H));
  num_cols = Math.ceil(canvas.height / (2 * H));
  for (let i = 0; i < num_rows * num_cols; i++) {
    grid.push([]);
  }
};

const createDamBreak = () => {
  for (let y = 60; y < canvas.height; y += H * 0.9)
    for (
      let x = canvas.width * (2 / 6);
      x <= canvas.width * (4 / 6);
      x += H * 0.9
    )
      if (particles.length < NUM_PARTICLES) {
        const jitter = Math.random() - 0.5;
        particles.push(new Particle(new Vector2(x + jitter, y)));
      }
};

const sort = () => {
  grid = [];
  initGrid();
  for (const p of particles) {
    //TODO: investigate why this is
    if (isNaN(p.position.x) || isNaN(p.position.y)) {
      return;
    }
    const row = Math.floor(p.position.x / (2 * H));
    const col = Math.floor(p.position.y / (2 * H));
    grid[row * num_cols + col].push(p);
  }
};

const adjacentCells = (i) => {
  const cells = [
    i,
    i + 1,
    i - 1,
    i - num_cols,
    i - num_cols + 1,
    i - num_cols - 1,
    i + num_cols,
    i + num_cols + 1,
    i + num_cols - 1,
  ]
    .filter((cell) => cell > 0 && cell < grid.length)
    .map((cell) => grid[cell]);
  return cells;
};

const computeDensityPressure = () => {
  for (let i = 0; i < grid.length; i++) {
    for (const pi of grid[i]) {
      pi.density = 0;
      for (const cell of adjacentCells(i)) {
        for (const pj of cell) {
          const rij = pj.position.sub(pi.position);
          const r2 = rij.len2();
          if (r2 < HSQ) {
            // this computation is symmetric
            pi.density += MASS * POLY6 * Math.pow(HSQ - r2, 3);
          }
        }
      }
      pi.pressure = GAS_CONST * (pi.density - REST_DENS);
    }
  }
};

const computeForces = () => {
  for (let i = 0; i < grid.length; i++) {
    for (const pi of grid[i]) {
      let fpress = new Vector2();
      let fvisc = new Vector2();
      let fst = new Vector2();
      for (const cell of adjacentCells(i))
        for (const pj of cell) {
          if (pi === pj) continue;
          const rij = pj.position.sub(pi.position);
          const r = rij.length();

          if (r < H && r > EPSILON) {
            const c1 =
              ((-MASS * (pi.pressure + pj.pressure)) / (2 * pj.density)) *
              SPIKY_GRAD *
              Math.pow(H - r, 2);
            fpress = fpress.add(rij.normalized().scale(c1));

            const c2 = ((VISC * MASS) / pj.density) * VISC_LAP * (H - r);
            fvisc = fvisc.add(pj.velocity.sub(pi.velocity).scale(c2));
          }
          if (r < 2 * H && r > EPSILON) {
            const c3 =
              -SIGMA *
              VISC_LAP *
              (MASS / pj.density) *
              SPIKY_GRAD *
              Math.pow(2 * H - r, 3);
            fst = fst.add(rij.normalized().scale(c3));
          }
        }

      const fgrav = G.scale(pi.density);
      // const fgrav = new Vector2();
      pi.force = fpress.add(fvisc).add(fgrav).add(fst);
    }
  }
};

const integrate = () => {
  for (let i = 0; i < particles.length; i++) {
    const p = particles[i];

    // forward Euler integration
    p.velocity = p.velocity.add(p.force.scale(DT / p.density));
    p.position = p.position.add(p.velocity.scale(DT));

    // enforce boundary conditions
    if (p.position.x < 0) {
      p.position.x = 0;
      p.velocity.x *= BOUND_DAMPING;
    }
    if (p.position.x > canvas.width - P_SIZE) {
      p.position.x = canvas.width - P_SIZE;
      p.velocity.x *= BOUND_DAMPING;
    }
    if (p.position.y < 0) {
      p.position.y = 0;
      p.velocity.y *= BOUND_DAMPING;
    }
    if (p.position.y > canvas.height - P_SIZE) {
      p.position.y = canvas.height - P_SIZE;
      p.velocity.y *= BOUND_DAMPING;
    }
  }
};

const drawParticles = () => {
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  particles.forEach((p) => {
    ctx.beginPath();
    ctx.arc(p.position.x, p.position.y, P_SIZE, 0, Math.PI * 2);
    ctx.fill();
    ctx.closePath();
  });
};

const flatten = () => {
  particles = [];
  for (const cell of grid) {
    for (const p of cell) particles.push(p);
  }
};

const main = () => {
  sort();
  computeDensityPressure();
  computeForces();
  flatten();
  integrate();
  drawParticles();
};
initGrid();
createDamBreak();
// main();
setInterval(main, 1);

document.addEventListener("keypress", ({ key }) => {
  if (key == " ") {
    const ps = [];
    for (let y = 60; y < canvas.height; y += H * 0.95)
      for (
        let x = canvas.width * (2 / 6);
        x <= canvas.width * (4 / 6);
        x += H * 0.9
      )
        if (ps.length < NUM_PARTICLES) {
          const jitter = Math.random() - 0.5;
          ps.push(new Particle(new Vector2(x + jitter, y)));
        }
    ps.forEach((p) => particles.push(p));
    sort();
  }
});
