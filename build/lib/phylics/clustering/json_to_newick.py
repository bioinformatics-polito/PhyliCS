def json_to_newick(json_obj):
	""" Return a newick string from a tree json format.
	Parameters:
	-----------
	json - A json object. See http://en.wikipedia.org/wiki/JSON
	Output:
	-------
	string - A string in newick format. See http://en.wikipedia.org/wiki/Newick_format
	Example:
	--------
	>>> json = {
		"name" : "F",
		"children": [
			{"name": "A", "branch_length": 0.1},
			{"name": "B", "branch_length": 0.2},
			{"name": "E","branch_length": 0.5,
			"children": [
				{"name": "C", "branch_length": 0.3},
				{"name": "D", "branch_length": 0.4}
				]
			}
		]
	}
	>>> print(json_to_newick(json))
	(A:0.1,B:0.2,(C:0.3,D:0.4)E:0.5)F;
	"""
	def _parse_json(json_obj):
		"""Return a json object in the format desribed below
		(------- daughter node ------,------- daughter node ------, ...,------- daughter node ------)-- root/parent -- 
		((children)name:branch_length,(children)name:branch_length, ...,(children)name:branch_length))name
		"""

		try:
			# Test is the key 'name' in the current level of the dictionary.
			newick = json_obj['name']
		except KeyError:
			# Catch no 'name' trees and start the newick string with empty qoutes
			newick = ''

		if 'edge_length' in json_obj:
			newick = newick + ':' + str(json_obj['edge_length'])

		# If there are 'children'	
		if 'children' in json_obj:
			# Initialise a list to contain the daughter info
			info = []
			# For each child, treat it as a new dictionary object
			for child in json_obj['children']:
				# parse the newick string straight into it the list
				info.append(_parse_json(child))

			# join all the daughter info together with a comma
			info = ','.join(info)

			# Concatenate all the children together at the start of the parent newick string		
			newick = '(' + info + ')' + newick

		return(newick)

	newick = _parse_json(json_obj) 

	return(newick)