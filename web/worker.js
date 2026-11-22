/* PathRx web worker: runs Python (Pyodide) off the page's main thread.

   Messages in:  {id, fn, args}        call pathrx_web.<fn>(*args)
                 {id, fn: "__init__", args: [config]}
   Messages out: {type: "status", text}
                 {type: "progress", label, loaded, total}
                 {id, ok: true, result} | {id, ok: false, error, traceback}
*/
let pyodide = null;
let web = null;

function post(msg) { self.postMessage(msg); }

const RETRY = new Set([403, 429, 500, 502, 503, 504]);  // NCBI answers 403 now and then under load
async function jsFetch(url, label) {
  let resp;
  for (let attempt = 0; ; attempt++) {
    try { resp = await fetch(url); } catch (e) { resp = null; if (attempt >= 3) throw e; }
    if (resp && (resp.ok || !RETRY.has(resp.status) || attempt >= 3)) break;
    post({type: "progress", label: `${label} (retrying)`, loaded: 0, total: 0});
    await new Promise((r) => setTimeout(r, 1500 * 2 ** attempt));
  }
  if (!resp.ok) throw new Error(`HTTP ${resp.status} for ${url}`);
  const total = Number(resp.headers.get("Content-Length")) || 0;
  if (!resp.body || !resp.body.getReader) {
    return new Uint8Array(await resp.arrayBuffer());
  }
  const reader = resp.body.getReader();
  const chunks = [];
  let loaded = 0, last = 0;
  for (;;) {
    const {done, value} = await reader.read();
    if (done) break;
    chunks.push(value);
    loaded += value.length;
    const now = Date.now();
    if (now - last > 150) { post({type: "progress", label, loaded, total}); last = now; }
  }
  post({type: "progress", label, loaded, total: total || loaded, done: true});
  const out = new Uint8Array(loaded);
  let off = 0;
  for (const c of chunks) { out.set(c, off); off += c.length; }
  return out;
}

function syncfs(populate) {
  return new Promise((resolve) => {
    try { pyodide.FS.syncfs(populate, () => resolve()); } catch (e) { resolve(); }
  });
}

async function init(config) {
  post({type: "status", text: "Loading Python (Pyodide)…"});
  importScripts(config.pyodide + "pyodide.js");
  pyodide = await loadPyodide({indexURL: config.pyodide});
  post({type: "status", text: "Loading NumPy, pandas, SciPy, statsmodels and matplotlib…"});
  await pyodide.loadPackage(["numpy", "pandas", "scipy", "statsmodels", "matplotlib"]);

  // downloads persist between visits in IndexedDB (falls back to memory)
  try {
    pyodide.FS.mkdir("/cache");
    pyodide.FS.mount(pyodide.FS.filesystems.IDBFS, {}, "/cache");
    await syncfs(true);
  } catch (e) { /* no IndexedDB: in-memory cache only */ }

  post({type: "status", text: "Loading PathRx…"});
  const zip = await (await fetch(config.bundle)).arrayBuffer();
  pyodide.unpackArchive(zip, "zip", {extractDir: "/lib"});
  pyodide.runPython("import sys; sys.path.insert(0, '/lib')");
  web = pyodide.pyimport("pathrx_web");
  web.js_fetch = jsFetch;
  const version = pyodide.runPython("import pathrx; pathrx.__version__ if hasattr(pathrx, '__version__') else ''");
  return JSON.stringify({pyodide: pyodide.version, pathrx: version});
}

async function clearCache() {
  pyodide.runPython(`
import shutil, os
for name in os.listdir('/cache'):
    p = os.path.join('/cache', name)
    shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
`);
  await syncfs(false);
  return "{}";
}

function cacheSize() {
  return pyodide.runPython(`
import os, json
files = [os.path.join('/cache', n) for n in os.listdir('/cache')] if os.path.isdir('/cache') else []
json.dumps({"files": len(files), "bytes": sum(os.path.getsize(f) for f in files if os.path.isfile(f))})
`);
}

self.onmessage = async (ev) => {
  const {id, fn, args = []} = ev.data;
  try {
    let result;
    if (fn === "__init__") result = await init(args[0]);
    else if (fn === "__clear_cache__") result = await clearCache();
    else if (fn === "__cache_size__") result = cacheSize();
    else {
      result = web[fn](...args);
      if (result && typeof result.then === "function") result = await result;
      if (result && typeof result.toJs === "function") result = result.toString();
      await syncfs(false);
    }
    post({id, ok: true, result});
  } catch (e) {
    const text = String(e && e.message || e);
    // Python errors arrive as a full traceback; the last line is the readable part
    const lines = text.trim().split("\n");
    post({id, ok: false, error: lines[lines.length - 1], traceback: text});
  }
};
