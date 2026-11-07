"""Fill the page template with the viz data and wrap it as a standalone HTML document.
Usage: python build_page.py template.html viz_data.json out.html"""
import sys
tpl, data, out = sys.argv[1:4]
html = open(tpl).read().replace("__DATA__", open(data).read().replace("</", "<\\/"))
if not html.lstrip().lower().startswith("<!doctype"):
    html = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover"></head><body>\n'
            + html + "\n</body></html>")
open(out, "w").write(html); print(out, round(len(html)/1e6, 2), "MB")
