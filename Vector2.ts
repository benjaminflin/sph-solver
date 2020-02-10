export default class Vector2 {
  x: number;
  y: number;
  constructor(x = 0, y = 0) {
    this.x = x;
    this.y = y;
  }

  add({ x, y }) {
    return new Vector2(this.x + x, this.y + y);
  }

  sub({ x, y }) {
    return new Vector2(this.x - x, this.y - y);
  }

  scale(s) {
    return new Vector2(this.x * s, this.y * s);
  }

  len2() {
    return this.x * this.x + this.y * this.y;
  }

  length() {
    return Math.sqrt(this.len2());
  }

  normalized() {
    const len = this.length();
    return new Vector2(this.x / len, this.y / len);
  }
}
