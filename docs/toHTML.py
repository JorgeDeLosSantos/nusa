# -*- coding: utf-8 -*-
import os
import subprocess

source = "docs/source"
build = "docs"
instr = "sphinx-build -b html "+source+" "+build
os.system(instr)
os.startfile(build+r"\index.html")

# instr = ["sphinx-build", "-b", "html", source, build]

# try:
#     # Ejecuta el comando
#     subprocess.run(instr, check=True)
#     # Abre el archivo index.html generado
#     os.startfile(os.path.join(build, "index.html"))
# except subprocess.CalledProcessError as e:
#     print(f"Error al ejecutar Sphinx: {e}")
