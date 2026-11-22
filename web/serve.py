"""Serve the PathRx web app locally, with KEGG included.

    python web/serve.py [--port 8000] [--pyodide DIR]

Builds web/dist/, serves it at http://localhost:8000/, and forwards
http://localhost:8000/kegg/<path> to https://rest.kegg.jp/<path>. Browsers block
direct requests to KEGG (it sends no CORS headers); going through this server
makes them same-origin, so the KEGG pathway sets, the KEGG DRUG route and network
proximity all run. Only rest.kegg.jp is forwarded, and only GET requests.

--pyodide DIR serves a local Pyodide distribution at /pyodide/ (for offline use);
open the page as http://localhost:8000/?pyodide=/pyodide/ then.
"""

from __future__ import annotations

import argparse
import functools
import http.server
import urllib.error
import urllib.request
from pathlib import Path

from build import WEB, build

KEGG = "https://rest.kegg.jp/"


class Handler(http.server.SimpleHTTPRequestHandler):
    pyodide_dir: Path | None = None

    def do_GET(self):
        if self.path.startswith("/kegg/"):
            return self._kegg(self.path[len("/kegg/"):])
        if self.pyodide_dir and self.path.startswith("/pyodide/"):
            rel = self.path.split("?", 1)[0][len("/pyodide/"):]
            f = (self.pyodide_dir / rel).resolve()
            if self.pyodide_dir in f.parents and f.is_file():
                self.send_response(200)
                ctype = {".js": "text/javascript", ".mjs": "text/javascript", ".wasm": "application/wasm",
                         ".json": "application/json"}.get(f.suffix, "application/octet-stream")
                self.send_header("Content-Type", ctype)
                self.send_header("Content-Length", str(f.stat().st_size))
                self.end_headers()
                self.wfile.write(f.read_bytes())
                return
            return self.send_error(404)
        return super().do_GET()

    def _kegg(self, rest: str):
        try:
            with urllib.request.urlopen(KEGG + rest, timeout=120) as r:
                body, status = r.read(), r.status
        except urllib.error.HTTPError as e:
            body, status = e.read(), e.code
        except urllib.error.URLError as e:
            return self.send_error(502, f"KEGG unreachable: {e.reason}")
        self.send_response(status)
        self.send_header("Content-Type", "text/plain; charset=utf-8")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--port", type=int, default=8000)
    ap.add_argument("--pyodide", type=Path, help="local Pyodide distribution to serve at /pyodide/")
    args = ap.parse_args()
    dist = build(WEB / "dist")
    Handler.pyodide_dir = args.pyodide.resolve() if args.pyodide else None
    handler = functools.partial(Handler, directory=str(dist))
    with http.server.ThreadingHTTPServer(("127.0.0.1", args.port), handler) as httpd:
        print(f"PathRx web app at http://localhost:{args.port}/  (KEGG forwarded at /kegg/)")
        httpd.serve_forever()


if __name__ == "__main__":
    main()
