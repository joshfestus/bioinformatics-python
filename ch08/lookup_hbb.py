# Ask Ensembl's REST API about one gene and print the JSON answer.
import json
import urllib.request

URL = "https://rest.ensembl.org/lookup/symbol/homo_sapiens/HBB"
request = urllib.request.Request(
    URL, headers={"Accept": "application/json"}
)

with urllib.request.urlopen(request) as response:
    answer = json.load(response)
print(json.dumps(answer, indent=2))
