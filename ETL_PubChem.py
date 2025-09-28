import logging
import pprint
import requests
from pymongo import MongoClient
from pymongo.errors import ConnectionFailure, PyMongoError
from rdkit import Chem
from rdkit.Chem import Descriptors

import logging

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def fetch_data(compound_name: str, base_url: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name"):
  """Fetch data """
  if not compound_name or not isinstance(compound_name, str):
    raise ValueError(f"Invalid compound name {compound_name}")
  url=f"{base_url}/{compound_name}/JSON"
  response = requests.get(url)
  logger.debug(f"Fetching from {url}")
  if response.status_code !=200:
    raise ValueError(f"Failed: {response.status_code}")
  
  data=response.json()
  pprint.pprint(data)
  if "PC_Compounds" not in data or not data["PC_Compounds"]:
    raise ValueError(f"No compound data found for: {compound_name}")
  compound = data["PC_Compounds"][0]
  #extract SMILEs
  props=compound.get("props")
  pprint.pprint(props)
  smiles=None
  for p in props:
    if p.get("urn", {}).get("label") == "SMILES" and p["urn"].get("name") == "Connectivity":
      smiles = p["value"]["sval"]
      break
  if not smiles:
    raise ValueError(f"No SMILES found for: {compound_name}")
  logger.info(f"Successfully fetched data for: {compound_name}")
  return {"name": compound_name, "smiles": smiles}

def compute_descriptors(smiles: str):
  """Compute molecular descriptors with RDKit"""
  if not smiles or not isinstance(smiles, str):
    raise ValueError(f"Pass SMILEs {smiles}")
  mol=Chem.MolFromSmiles(smiles)
  if mol is None:
    raise ValueError("Failed to transform SMILES into molecules ")
  return {
    "MolWt": Descriptors.MolWt(mol),
    "LogP": Descriptors.MolLogP(mol),
    "NumHDonors": Descriptors.NumHDonors(mol),
    "NumHAcceptors": Descriptors.NumHAcceptors(mol),
  }
  
def save_to_mongo(document: dict, db_name="science_db", collection_name="compounds"):
  """Save to MongoDB"""
  try:
    client = MongoClient("mongodb://localhost:27017")
    db = client[db_name]
    collection = db[collection_name]
    result = collection.insert_one(document)
    return result.inserted_id
  except (PyMongoError, ConnectionFailure) as e:
    logger.error(f"MongoDB error: {e}")
    raise
  finally:
    client.close()
    

if __name__=="__main__":
  try:
    #Extract
    compound="aspirin"
    pubchem_data =fetch_data(compound)
    print("Fetched from PubChem:", pubchem_data)
  
    # Transform
    descriptors = compute_descriptors(pubchem_data["smiles"])
    print("Descriptors:", descriptors)
  
    # Load
    doc = {**pubchem_data, "descriptors": descriptors}
    mongo_id = save_to_mongo(doc)
    print(f"Document saved to MongoDB с _id={mongo_id}")
  except ValueError as ve:
    logger.error(f"Validation error: {ve}")
  except requests.RequestException as re:
    logger.error(f"Network error: {re}")
  except PyMongoError as pe:
    logger.error(f"MongoDB error: {pe}")
  except Exception as e:
    logger.critical(f"Unexpected error: {e}")
    