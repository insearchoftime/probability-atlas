import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import vm from 'node:vm';
import {fileURLToPath} from 'node:url';

const root = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '..');
const site = path.join(root, '_site');
const html = fs.readFileSync(path.join(site, 'index.html'), 'utf8');
let count = 0;
for (const match of html.matchAll(/<script\b([^>]*)>([\s\S]*?)<\/script>/gi)) {
  const source = /\bsrc=["']([^"']+)["']/.exec(match[1]);
  const name = source ? source[1] : `inline-script-${++count}`;
  const code = source ? fs.readFileSync(path.join(site, name), 'utf8') : match[2];
  new vm.Script(code, {filename: name});
}
for (const match of html.matchAll(/<link\b[^>]*href=["']([^"']+\.css)["']/gi)) {
  if (!match[1].startsWith('http')) assert.ok(fs.existsSync(path.join(site, match[1])), match[1]);
}
assert.ok(html.includes("'t-statistic','z-statistic'"), 'Both new cards must be on the visible atlas route');
assert.ok(html.includes('atlas-editorial.js'), 'Editorial changes must be loaded');
assert.ok(html.includes('atlas-statistics.js'), 'Statistics must be loaded');
console.log(`All ${count} inline scripts and companion scripts parse; publication assets and route are present.`);
