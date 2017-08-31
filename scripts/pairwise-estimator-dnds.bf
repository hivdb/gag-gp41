RequireVersion ("2.3.2");
LoadFunctionLibrary ("libv3/all-terms.bf");
LoadFunctionLibrary ("libv3/tasks/alignments.bf");
LoadFunctionLibrary ("libv3/tasks/genetic_code.bf");
LoadFunctionLibrary ("libv3/tasks/estimators.bf");
LoadFunctionLibrary ("libv3/UtilityFunctions.bf");
LoadFunctionLibrary ("libv3/convenience/regexp.bf");

LoadFunctionLibrary ("libv3/models/codon/MG_REV.bf");

pairwise.alignment  = alignments.ReadCodonDataSet ("pairwise.dataset");

alignments.LoadCodonDataFile ("pairwise.dataset", "pairwise.filter", pairwise.alignment);

pairwise.names      = alignments.GetSequenceNames ("pairwise.dataset");

pairwise.groups     = {};

utility.ForEach (pairwise.names, "_name_", "
        key = pairwise.name2key (_name_);
        utility.EnsureKey (pairwise.groups, key);
        pairwise.groups[key] + _name_;

    ");

lfunction pairwise.name2key (name) {
    parts = regexp.Split (name, '_');
    // remove last
    parts - (Abs (parts) - 1);
    return Join ("_", parts);
}


alignments.GetSequenceByName ("pairwise.dataset", None);

pairwise.filter_string = io.PromptUserForString ("Specify an optional, 0-based, filter string for codons (e.g. 0-100); enter 'all' to include all codons");

if (pairwise.filter_string == "all") {
    pairwise.filter_string = "";
}

console.log ("\nPair, dN/dS");
pairwise.results = LAST_FILE_PATH;

pairwise.filters = {};
pairwise.trees = {};
pairwise.index_to_name = {};

utility.ForEachPair (pairwise.groups, "_key_", "_sequences_",
'
    this_pair = pairwise.handle_pair (_key_, _sequences_);
    pairwise.filters + this_pair ["filter"];
    pairwise.trees   + this_pair ["tree"];
    pairwise.index_to_name [Abs (pairwise.trees) - 1] = _key_;
');


pairwise.mles = estimators.FitMGREV (pairwise.filters, pairwise.trees, (^"pairwise.alignment")["code"], {"model-type" : ^"terms.local"}, None);

pairwise.results = {};

utility.ForEachPair (pairwise.mles[terms.branch_length], "_key_", "_value_",
'
	dn_ds = utility.Map (_value_, "_mles_", "(_mles_[terms.parameters.nonsynonymous_rate])[terms.fit.MLE]/(_mles_[terms.parameters.synonymous_rate])[terms.fit.MLE]");
	console.log (pairwise.index_to_name[_key_] + ", " + (utility.Values(dn_ds))[0]);
	pairwise.results + (utility.Values(dn_ds))[0];
');


//console.log ("c (" + Join (",", pairwise.results) + ")");



lfunction pairwise.handle_pair (key, data) {
    assert (Abs (data) == 2, "Expected exactly two sequences per sequence group. Had " + Abs (data) + " for `key`");

    filter = {data[0] : 1, data[1]: 1};

    utility.ExecuteInGlobalNamespace ('
    	DataSetFilter `&this_pair`.`key` = CreateFilter (^"pairwise.filter", 3, "`^"pairwise.filter_string"`", `&filter`[(^"pairwise.names")[speciesIndex]], (^"pairwise.alignment")["stop"]);
    ');

	tree_string = "(" + key[0] + "," + key[1] + ")";

	return {"filter" : "`&this_pair`.`key`", "tree": {"string" : tree_string}};


}
