import sys
tpl, data, out = sys.argv[1:4]
html = open(tpl).read().replace("__DATA__", open(data).read().replace("</", "<\\/"))
open(out, "w").write(html); print(out, round(len(html)/1e6, 2), "MB")
