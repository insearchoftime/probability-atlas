'use strict';
const assert = require('node:assert/strict');
require('../atlas-statistics.js');
const {math, codeFor} = global.AtlasStatistics;
const close = (actual, expected, tolerance = 2e-10) => assert.ok(Math.abs(actual - expected) <= tolerance, `${actual} != ${expected}`);

// Independent published distribution values; these catch wrong tails, df and SE.
close(math.normalTwoSidedP(2), 0.04550026389635839);
close(math.tTwoSidedP(2, 4), 0.1161165235168155);
close(math.tTwoSidedP(2, 15), 0.0639450072847202);
close(math.tTwoSidedP(1e-6, 9999), 0.9999992021353881, 1e-10);
close(math.chiP(4, 1), 0.04550026389635857);
close(math.inverseTail(.05, x => math.tTwoSidedP(x, 4)), 2.7764451051977987);
close(math.inverseTail(.05, x => math.chiP(x, 5)), 11.070497693516351);

const p = {mean:102, mu:100, sd:4, n:16, alpha:.05};
const z = math.meanModel('z-statistic', p), t = math.meanModel('t-statistic', p);
close(z.se,1); close(z.statistic,2); assert.equal(z.reject,true);
assert.equal(t.df,15); assert.equal(t.reject,false); assert.ok(t.lower < p.mu);
assert.ok(z.lower > p.mu); assert.ok(t.p > z.p);

// Same relative departure gives a larger chi-square as independent N grows.
const a = math.chiModel([9,15,10,16,7,3]);
const b = math.chiModel([90,150,100,160,70,30]);
close(a.statistic,12); close(b.statistic,120); assert.ok(b.p < a.p);
close(math.chiModel([10,10,10,10,10,10]).p,1);
assert.equal(math.chiModel([1,1,1,1,1,1]).adequate,false);
assert.equal(Number.isNaN(math.chiModel([0,0,0,0,0,0]).p),true);

for (const df of [1,2,4,15,100,9999]) {
  let previous = 1;
  for (const x of [0,.001,.1,1,2,5,10,100,10000]) {
    const p = math.tTwoSidedP(x,df);
    assert.ok(p >= 0 && p <= previous); previous = p;
    close(p,math.tTwoSidedP(-x,df));
  }
  for(const alpha of [.1,.05,.01]){
    const c = math.inverseTail(alpha,x=>math.tTwoSidedP(x,df));
    close(math.tTwoSidedP(c,df),alpha,1e-11);
  }
}
for(const id of ['z-statistic','t-statistic','pearson']) {
  const code = codeFor(id);
  assert.ok(code.code.includes('Текущие параметры карточки'));
  assert.ok(code.filename.endsWith('.py'));
}
console.log('Atlas statistics: reference values, interval/test agreement, chi-square scaling, tail monotonicity and symmetry passed.');
