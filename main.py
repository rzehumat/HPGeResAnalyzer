import rptParser
import getIg
import getEpsilon

from sys import argv

rptParser.parse_RPT(argv[1])
getIg.append_Igamma_dir("parsed_reports/")
getEpsilon.add_epsilon_dir("with_Ig/")