# Copyright 2013 by Kamil Koziara. All rights reserved.               
# This code is part of the Biopython distribution and governed by its    
# license.  Please see the LICENSE file that should have been included   
# as part of this package.

"""
Module containing classes for finding functional gene enrichments using
ontologies.
"""
from __future__ import print_function
from Bio._py3k import range, filter


import collections
import random, math, bisect
from . import Stats
from . import IdResolver

class EnrichmentEntry(object):
    """
    Represents one result returned by TermForTermEnrichmentFinder.
    """
    def __init__(self, term_id, term_name, p_value):
        self.id = term_id
        self.name = term_name
        self.p_value = p_value
        self.corrections = {}
        self.attrs = {}
    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __repr__(self):
        return 'EnrichmentEntry(id = {0}, name = {2}, p_value = {1})'.format(self.id, self.p_value, self.name)

    def __str__(self):
        return """ID : {0}
name : {1}
p-value : {2}
corrected p-values: {3}""".format(self.id, self.name, self.p_value, self.corrections)

class Enrichment(object):
    """
    Contains all results found by TermForTermEnrichmentFinder
    """
    
    def __init__(self, method, entries, warnings, corrections = []):
        self.entries = entries
        self.warnings = warnings
        self.method = method
        self.corrections = corrections
    
    def filter(self, filter_fun):
        return Enrichment(self.method, 
            list(filter(filter_fun, self.entries)),
            list(self.warnings), list(self.corrections))
        
    def filter_p_val(self, p_val):
        """
        Returns enrichment entries with specified siginificance level.
        """
        return self.filter(lambda x: x.p_value <= p_val)
    
    def __eq__(self, other):
        if isinstance(other, self.__class__):
            return self.__dict__ == other.__dict__
        else:
            return False

    def __ne__(self, other):
        return not self.__eq__(other)
    
    def __repr__(self):
        return "Enrichment(method = {2}, entries_num = {0} , warnings_num = {1})".format(len(self.entries), len(self.warnings), self.method)

    def __str__(self):
        return "Enrichment found using {2} method: {0} entries, {1} warnings.".format(len(self.entries),
                                                               len(self.warnings), self.method)
class EnrichmentFinder(object):
    
    def find_enrichment(self, genes, **args):
        raise NotImplementedError(("find_enrichment not implemented. ", 
                                   "This is just a base class."))    
    
class BaseEnrichmentFinder(EnrichmentFinder):
    
    def __init__(self, annotations, ontology_graph, resolver_generator):
        """
        Initialize base for Enrichment Finder.
        """
        
        self.ontology_graph = ontology_graph
        self.annotations = annotations
        self.resolver = resolver_generator(iter(self.annotations.values()))
    

    def _find_terms_associations(self, gene_list):
        terms_assocs = collections.defaultdict(set)
        for gene in gene_list:
            if gene in self.annotations:
                enriched_terms = set()
                for term in self.annotations[gene].associations:
                    node = self.ontology_graph.get_node(term.term_id)
                    if node != None:
                        nid = node.data.id # because enrichments may use synonyms instead of ids
                        enriched_terms.add(nid)
                        enriched_terms |= self.ontology_graph.get_ancestors(nid)
                for t in enriched_terms:
                    terms_assocs[t].add(gene)
        return terms_assocs

    def _resolve_ids(self, gene_list, warnings):
        resolved_list = []
        for x in gene_list:
            rx = self.resolver.resolve(x)
            if x != rx:
                warnings.append("Unknown id: '{0}' was resolved to: '{1}'".format(x, rx))
            resolved_list.append(rx)
        return resolved_list
    
    @staticmethod
    def _calculate_corrections(result, corrections):
        """
        Calculates corrections.
        """
        
        if len(corrections) > 0:
            pvals = [x.p_value for x in result]
            corr_pvals = []
            for c_id in corrections:
                cfun = Stats.corrections[c_id]
                corr_pvals.append((c_id, cfun(pvals)))
            for i in range(len(result)):
                result[i].corrections = dict([(c_id, pv[i]) for c_id, pv in corr_pvals])
    
class TermForTermEnrichmentFinder(BaseEnrichmentFinder):
    """
    Utility for finding enriched group of terms given list of genes connected
    to the tested phenotype, ontology graph and associations to this graph.
    
    TermForTermEnrichmentFinder could be used for finding enrichment in many cases of
    entities but in this example let's focus on finding gene enrichment in
    gene ontology graph.

    To initialize the finder you need at least ontology graph and associations
    from genes to it's nodes. Typically you would read them from files
    using readers from Bio.Ontology.IO module.
    
    >>> import Bio.Ontology.IO as OntoIO
    >>> go_graph = OntoIO.read("Ontology/go_test.obo", "obo")
    >>> assocs = OntoIO.read("Ontology/ga_test.fb", "gaf")
    
    Now is the time to create TermForTermEnrichmentFinder. Besides the arguments
    mentioned before you could specify an id resolver, which basically tries
    to find synonyms for gene ids that finder can understand, and population
    of genes as reference. If none of this is specified default resolver
    is used and all genes from association are used as the population.
    
    >>> from Bio.Ontology import TermForTermEnrichmentFinder
    >>> ef = TermForTermEnrichmentFinder(assocs, go_graph)
    
    To run finder you just need to call find_enrichment method with list
    of genes as an argument:
    
    >>> genes_to_study = ['FBgn0070057', '18-wheeler', 'FBgn0043467']
    
    Additionally you can add a list of corrections which should be used
    (by default it's empty):
    
    >>> corrections = ['bonferroni', 'bh_fdr']
    >>> result = ef.find_enrichment(genes_to_study, corrections)
    
    The result contains a list of EnrichmentEntry instances - each for a node in graph
    which is enriched - and a list of warnings.
    
    >>> print(result)
    Enrichment found using term_for_term method: 64 entries, 1 warnings.
    
    Notice there is one warning. Let's print the list:
    >>> print(result.warnings)
    ["Unknown id: '18-wheeler' was resolved to: 'FBgn0004364'"]
    
    Earlier when specifying list of genes we used non-standardized gene id:
    '18-wheeler'. Thanks to default resolver (FirstOneResolver)
    of TermForTermEnrichmentFinder standard id was inferred.
    
    >>> print(result.entries[0])
    ID : GO:0044707
    name : single-multicellular organism process
    p-value : 1.0
    corrected p-values: {'bh_fdr': 1.0, 'bonferroni': 1.0}
    
    """
    
    def __init__(self, annotations, ontology_graph, population = None,
                 resolver_generator = IdResolver.FirstOneResolver):
        """
        Initialize Enrichment Finder.
        
        Parameters
        ----------
        annotations - iterable containing annotations
        ontology_graph - graph with ontology
        population - population used as reference to study sample
        resolver_generator - constructor of resolver used to disambiguate ids
        """
        super(TermForTermEnrichmentFinder, self).__init__(annotations, ontology_graph, resolver_generator)

        if population != None:
            self.population = population
        else:
            self.population = list(self.annotations.keys())
        self.terms_to_population_genes = self._find_terms_associations(self.population)

        
    def find_enrichment(self, gene_list, corrections = []):
        """
        Finds enrichment of specified group of genes.
        
        Parameters
        ----------
        gene_list - list of genes to study
        corrections - list of corrections that should be applied to result.
            Possible values are:
                o "bonferroni" - Bonferroni correction,
                o "bh_fdr" - Benjamin-Hochberg FDR correction.
        """
        
        result = []
        warnings = []
        
        resolved_list = self._resolve_ids(gene_list, warnings)
            
        terms_to_study_genes = self._find_terms_associations(resolved_list)
        
        
        study_size = len(resolved_list)

        population_size = len(self.population)
        # Calculate enrichment for every term given the counts (term for term)
        for term, study_set in terms_to_study_genes.items():
            study_hits = len(study_set)
            population_hits = len(self.terms_to_population_genes[term])
            pval = Stats.hypergeometric_test(study_hits, study_size,
                                             population_hits, population_size)

            entry = EnrichmentEntry(term, self.ontology_graph.get_term(term).name, pval)

            result.append(entry)
        # Calculate chosen corrections
        BaseEnrichmentFinder._calculate_corrections(result, corrections)
        
        # check for warnings
        if len(self.ontology_graph.cycles) > 0:
            warnings.append("Graph contains cycles: " + str(self.ontology_graph.cycles))
        return Enrichment("term_for_term", result, warnings, corrections)


class ParentChildEnrichmentFinder(BaseEnrichmentFinder):
    """
    Finder implementing different method for finding enrichment called
    parent-child.
    
    Parent-child method takes relationships between the nodes into account.
    You can use it exactly like the term for term method:
    
    >>> import Bio.Ontology.IO as OntoIO
    >>> go_graph = OntoIO.read("Ontology/go_test.obo", "obo")
    >>> assocs = OntoIO.read("Ontology/ga_test.fb", "gaf")

    >>> from Bio.Ontology import ParentChildEnrichmentFinder
    >>> ef = ParentChildEnrichmentFinder(assocs, go_graph)
    >>> genes_to_study = ['FBgn0070057', '18-wheeler', 'FBgn0043467']
    >>> corrections = ['bonferroni', 'bh_fdr']
    
    >>> result = ef.find_enrichment(genes_to_study, corrections, "intersection")

    
    >>> print(result)
    Enrichment found using parent_child_intersection method: 61 entries, 1 warnings.
    >>> print(result.entries[0])
    ID : GO:0044707
    name : single-multicellular organism process
    p-value : 1.0
    corrected p-values: {'bh_fdr': 1.0, 'bonferroni': 1.0}
    
    """
    def __init__(self, annotations, ontology_graph, population = None,
                 resolver_generator = IdResolver.FirstOneResolver):
        """
        Initialize Enrichment Finder.
        
        Parameters
        ----------
        annotations - iterable containing annotations
        ontology_graph - graph with ontology
        population - population used as reference to study sample
        resolver_generator - constructor of resolver used to disambiguate ids
        """
        super(ParentChildEnrichmentFinder, self).__init__(annotations, ontology_graph, resolver_generator)

        if population != None:
            self.population = population
        else:
            self.population = list(self.annotations.keys())
        self.terms_to_population_genes = self._find_terms_associations(self.population)

    
    def _count_op_items(self, set_list, set_op):
        return len(set_op(*set_list)) if len(set_list) > 0 else 0
        
    def find_enrichment(self, gene_list, corrections = [], method = "union"):
        """
        Finds enrichment of specified group of genes. Method takes
        the parent-child relationship into account when computing p-value.
        Reference: http://bioinformatics.oxfordjournals.org/content/23/22/3024.long
        
        Parameters
        ----------
        gene_list - list of genes to study
        corrections - list of corrections that should be applied to result.
            Possible values are:
                o "bonferroni" - Bonferroni correction,
                o "bh_fdr" - Benjamin-Hochberg FDR correction.
        method - method of computing the p-values
            Possible values are:
                o "union",
                o "intersection"
        """
        
        result = []
        warnings = []
        
        resolved_list = self._resolve_ids(gene_list, warnings)
            
        terms_to_study_genes = self._find_terms_associations(resolved_list)

        for term, study_set in terms_to_study_genes.items():
            study_hits = len(study_set)
            population_hits = len(self.terms_to_population_genes[term])
            study_set_list = []
            pop_set_list = []
            # calculate sets of genes annotated to parents
            for parent in self.ontology_graph.get_parents(term):
                study_set_list.append(terms_to_study_genes[parent])
                pop_set_list.append(self.terms_to_population_genes[parent])
            
            if method == "union":
                set_op = set.union
            elif method == "intersection":
                set_op = set.intersection
            else:
                raise ValueError("{0} is not correct method type.".format(method))
            
            parents_study_size = self._count_op_items(study_set_list, set_op)
            population_parents_size = self._count_op_items(pop_set_list, set_op)
            
            if study_hits <= parents_study_size and population_hits <= population_parents_size:
                pval = Stats.hypergeometric_test(study_hits, parents_study_size,
                                             population_hits, population_parents_size)

                entry = EnrichmentEntry(term, self.ontology_graph.get_term(term).name, pval)
                entry.attrs = {"study_hits" : study_hits, "parents_study_size": parents_study_size,
                            "population_hits" : population_hits, "population_parents_size" : population_parents_size}
                result.append(entry)
        
        # Calculate chosen corrections
        BaseEnrichmentFinder._calculate_corrections(result, corrections)
        
        # check for warnings
        if len(self.ontology_graph.cycles) > 0:
            warnings.append("Graph contains cycles: " + str(self.ontology_graph.cycles))
        return Enrichment("parent_child_" + method, result, warnings, corrections)


class GseaEnrichmentFinder(BaseEnrichmentFinder):
    """
    Utility for finding enriched group of terms given list of genes ranked
    by correlation with given phenotype, ontology graph and associations to this graph.

    The usage of GseaEnrichmentFinder is very similar to TermForTermEnrichmentFinder.
    You have to initialize the finder with ontology graph and list of
    associations:
    
    >>> import Bio.Ontology.IO as OntoIO
    >>> go_graph = OntoIO.read("Ontology/go_test.obo", "obo")
    >>> assocs = OntoIO.read("Ontology/ga_test.fb", "gaf")
    >>> from Bio.Ontology import GseaEnrichmentFinder
    >>> ef = GseaEnrichmentFinder(assocs, go_graph)
    
    To run finder you just need to call find_enrichment
    method with gene rank:
    
    >>> genes_rank = [('FBgn0070057', 0.8), ('18-wheeler', 0.6), ('FBgn0043467', 0.2), ('FBgn0004222', -0.5)]
    
    Additionally you can specify permutations number (more is better, but slower).
    
    >>> result = ef.find_enrichment(genes_rank, perms_no = 1000)
    
    The result contains a list of EnrichmentEntry instances - each for a node in graph
    which is enriched - and a list of warnings.
    
    >>> print(result)
    Enrichment found using GSEA method: 12 entries, 1 warnings.
    
    """
    
    def __init__(self, annotations, ontology_graph,
                 resolver_generator = IdResolver.FirstOneResolver):
        """
        Initialize Enrichment Finder.
        
        Parameters
        ----------
        annotations - iterable containing annotations
        ontology_graph - graph with ontology
        resolver_generator - constructor of resolver used to disambiguate ids
        """
        
        super(GseaEnrichmentFinder, self).__init__(annotations, ontology_graph,
                                                     resolver_generator)
    def _get_perms(self, gene_list, perms_no):
        perms = []
        permutation = list(gene_list)
        for _ in range(perms_no):
            random.shuffle(permutation)
            perms.append(list(permutation))
        return perms
    
    def find_enrichment(self, gene_rank, perms_no = 1000,
                        min_set_rank_intersection = 2,  corr_power = 1.):
        """
        Finds enrichment using GSEA method.
        
        Reference: http://www.pnas.org/content/102/43/15545.full
        
        Parameters
        ----------
        - gene_rank
        - perms_no - number of permutations used to compute p-value
        - min_set_rank_intersection - minimal number of genes common to
          the set and rank to take the set into account
        - corr_power - weight of correlation when computing enrichment score
        
        """
        
        sorted_gene_rank = sorted(gene_rank, key = lambda x: x[1], reverse = True)
        gene_corr = []
        gene_list = []
        
        for pos in sorted_gene_rank:
            gene_list.append(pos[0])
            gene_corr.append(pos[1])
        
        warnings = []
        result = []
        resolved_list = self._resolve_ids(gene_list, warnings)
        enriched_terms = self._find_enriched_terms(resolved_list, min_set_rank_intersection)
        perms = self._get_perms(resolved_list, perms_no)
        
        # Computing both: uncorrected and FDR corrected p-value
        all_pos_nes_perm = []
        all_neg_nes_perm = []
        all_nes = []
        
        for term, gene_set in enriched_terms.items():
            
            orig_es, orig_plot = Stats.kolmogorov_smirnov_rank_test(gene_set, resolved_list, gene_corr, corr_power)
            
            pcount = 0
            pos_nes_perm = []
            neg_nes_perm = []
            
            for perm in perms:
                perm_es, _ = Stats.kolmogorov_smirnov_rank_test(gene_set, perm, gene_corr, corr_power)
                if orig_es < 0:
                    pcount += int(orig_es >= perm_es)
                else:
                    pcount += int(orig_es <= perm_es)
                
                # For FDR
                if perm_es < 0:
                    neg_nes_perm.append(perm_es)
                else:
                    pos_nes_perm.append(perm_es)
            
            # Computing p-value
            pval = pcount / float(perms_no)
            entry = EnrichmentEntry(term, self.ontology_graph.get_term(term).name,
                                    pval)
            entry.attrs["score"] = orig_es
            entry.attrs["plot"] = orig_plot
            result.append(entry)
            
            # Now part of the FDR computation
            neg_mean = Stats.mean(neg_nes_perm)
            pos_mean = Stats.mean(pos_nes_perm)
            
            if neg_mean != .0:
                neg_nes_perm = [- x / neg_mean for x in neg_nes_perm]
                if orig_es < 0:
                    orig_es = - orig_es / neg_mean
            if pos_mean != .0:
                pos_nes_perm = [x / pos_mean for x in pos_nes_perm]
                if orig_es > 0:
                    orig_es = orig_es / pos_mean
            
            all_nes.append(orig_es)
            all_neg_nes_perm += neg_nes_perm
            all_pos_nes_perm += pos_nes_perm
            
        # Last stage of FDR computation
        all_neg_nes = []
        all_pos_nes = []

        for nes in all_nes:
            if nes < 0:
                all_neg_nes.append(nes)
            else:
                all_pos_nes.append(nes)
        
        all_pos_nes_perm.sort()
        all_neg_nes_perm.sort()
        all_pos_nes.sort()
        all_neg_nes.sort()
        
        for i in range(len(all_nes)):
            nes = all_nes[i]
            if nes < 0:
                a = bisect.bisect_right(all_neg_nes_perm, nes) / float(len(all_neg_nes_perm))
                b = (bisect.bisect_right(all_neg_nes, nes) - 1) / float(len(all_neg_nes))
            else:
                a = 1 - bisect.bisect_left(all_pos_nes_perm, nes) / float(len(all_pos_nes_perm))
                b = 1 - (bisect.bisect_left(all_pos_nes, nes) + 1) / float(len(all_pos_nes))
            fdr = a / b if b != .0 else 1.
            if fdr > 1.:
                fdr = 1.
            result[i].corrections['fdr'] = fdr
            
        if len(self.ontology_graph.cycles) > 0:
            warnings.append("Graph contains cycles: " + str(self.ontology_graph.cycles))
        return Enrichment("GSEA", result, warnings, ['fdr'])
        
    def _find_enriched_terms(self, gene_list, min_set_size):
        enriched_terms = {}
        term_assocs = self._find_terms_associations(gene_list)
        for k, v in term_assocs.items():
            if len(v) >= min_set_size:
                enriched_terms[k] = v
        return enriched_terms

class RankedParentChildEnrichmentFinder(BaseEnrichmentFinder):
    """
    Utility for finding enriched group of terms given list of genes ranked
    by correlation with given phenotype, ontology graph and associations to this graph.

    This one uses novell method for finding enrichments - ranked parent-child:
    
    >>> import Bio.Ontology.IO as OntoIO
    >>> go_graph = OntoIO.read("Ontology/go_test.obo", "obo")
    >>> assocs = OntoIO.read("Ontology/ga_test.fb", "gaf")
    >>> from Bio.Ontology import RankedParentChildEnrichmentFinder
    >>> ef = RankedParentChildEnrichmentFinder(assocs, go_graph)
    
    >>> genes_rank = [('FBgn0070057', 0.8), ('18-wheeler', 0.6), ('FBgn0043467', 0.2), ('FBgn0004222', -0.5)]

    >>> result = ef.find_enrichment(genes_rank)
    >>> print(result)
    Enrichment found using ranked parent-child method: 61 entries, 1 warnings.
    
    """
    def __init__(self, annotations, ontology_graph,
                 resolver_generator = IdResolver.FirstOneResolver):
        """
        Initialize Enrichment Finder.
        
        Parameters
        ----------
        annotations - dictionary containing annotations
        ontology_graph - graph with ontology
        resolver_generator - constructor of resolver used to disambiguate ids
        """
        
        super(RankedParentChildEnrichmentFinder, self).__init__(annotations, ontology_graph,
                                                     resolver_generator)
    
    def _get_half_results(self, resolved_list, ef, method, warnings):
        list_slice = []
        results = collections.defaultdict(list)

        for gene in resolved_list:
            list_slice.append(gene)
            slice_res = ef.find_enrichment(list_slice, method = method)
            warnings += slice_res.warnings
            for e in slice_res.entries:
                results[e.id].append(e.p_value)
        return results
    
    def find_enrichment(self, gene_rank, side = "+", corrections = [],
                                     rank_as_population = False, method = "union"):
        """
        Finds enrichment by applying parent-child analysis to list slices.
        
        Parameters
        ----------
        - gene_rank
        - side - states which side of the rank (ordered by correlation) we are interested in
          o "+" - highest correlation
          o "-" - lowest correlation
          o "+/-" - both
        - corrections - corrections that shuld be applied
        - rank_as_population - if set to True only the genes in the rank will be set as
            the population,
        - method - method of parent-child to use
          o "union"
          o "intersection"
          
        """
        
        warnings = []
        
        sorted_gene_rank = sorted(gene_rank, key = lambda x: x[1], reverse = True)
        gene_list = [g for g, _ in sorted_gene_rank]
        resolved_list = self._resolve_ids(gene_list, warnings)
        
        
        if rank_as_population:
            ef = ParentChildEnrichmentFinder(self.annotations, self.ontology_graph,
                                  resolved_list, IdResolver.Resolver)
        else:
            ef = ParentChildEnrichmentFinder(self.annotations, self.ontology_graph,
                                  resolver_generator = IdResolver.Resolver)
        
        if side == "-":
            all_results = self._get_half_results(resolved_list[::-1], ef, method, warnings)
        elif side == "+":
            all_results = self._get_half_results(resolved_list, ef, method, warnings)
        elif side == "+/-":
            minus_results = self._get_half_results(resolved_list[::-1], ef, method, warnings)
            all_results = self._get_half_results(resolved_list, ef, method, warnings)
            for k, v in minus_results.items():
                all_results[k] += v
        else:
            raise ValueError('"{0}" is not correct side specification.'.format(side))
        
        result = []
        for oid, p_vals in all_results.items():
            min_pval = 1.0
            plot = []
            for pv in p_vals:
                if pv < min_pval:
                    min_pval = pv
                plot.append(1 - pv)
            entry = EnrichmentEntry(oid, self.ontology_graph.get_term(oid).name,
                                          min_pval)
            entry.attrs["plot"] = plot
            entry.attrs["score"] = 1 - min_pval
            
            result.append(entry)
            
        BaseEnrichmentFinder._calculate_corrections(result, corrections)
        
        if len(self.ontology_graph.cycles) > 0:
            warnings.append("Graph contains cycles: " + str(self.ontology_graph.cycles))
            
        return Enrichment("ranked parent-child", result, warnings, corrections)

if __name__ == "__main__":
    from Bio._utils import run_doctest
    run_doctest()

