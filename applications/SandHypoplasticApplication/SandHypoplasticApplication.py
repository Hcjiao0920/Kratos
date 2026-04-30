from KratosMultiphysics import _ImportApplication
from KratosSandHypoplasticApplication import *

application = KratosSandHypoplasticApplication()
application_name = "KratosSandHypoplasticApplication"

_ImportApplication(application, application_name)