import rptParser
import getIg
import getEpsilon
import os

from sys import argv
from shutil import move
from pathlib import Path

folder = argv[1]
geom = argv[2]

rptParser.parse_RPT(folder)
getIg.append_Igamma_dir("parsed_reports/")
getEpsilon.add_epsilon_dir("with_Ig", geom)

Path("out").mkdir(parents=True, exist_ok=True)

source_dir = "with_eps/"
files = os.listdir(source_dir)

for file_name in files:
    move(source_dir + file_name, "out")