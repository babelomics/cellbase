package org.opencb.cellbase.lib.impl;

import org.bson.Document;
import org.bson.conversions.Bson;
import org.hamcrest.CoreMatchers;
import org.junit.Before;
import org.junit.Test;
import org.opencb.biodata.models.core.Region;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.cellbase.core.api.ClinicalDBAdaptor;
import org.opencb.cellbase.core.loader.LoadRunner;
import org.opencb.cellbase.lib.GenericMongoDBAdaptorTest;
import org.opencb.commons.datastore.core.Query;
import org.opencb.commons.datastore.core.QueryOptions;
import org.opencb.commons.datastore.core.QueryResult;

import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.net.URISyntaxException;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.ExecutionException;
import java.util.regex.Pattern;
import java.util.stream.Collectors;

import static junit.framework.TestCase.assertNotNull;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertThat;
import static org.junit.Assert.assertTrue;

/**
 * Created by fjlopez on 24/03/17.
 */
public class ClinicalMongoDBAdaptorTest extends GenericMongoDBAdaptorTest {

    public ClinicalMongoDBAdaptorTest() throws IOException {
    }

    @Before
    public void setUp() throws Exception {
        clearDB(GRCH37_DBNAME);
        Path path = Paths.get(getClass()
                .getResource("/clinical_variants.full.test.json.gz").toURI());
        loadRunner.load(path, "clinical_variants");
        Path pathDiseases = Paths.get(getClass()
                .getResource("/diseases.test.json.gz").toURI());
        loadRunner.load(pathDiseases, "diseases");
    }

    @Test
    public void nativeGet() throws Exception {
        ClinicalDBAdaptor clinicalDBAdaptor = dbAdaptorFactory.getClinicalDBAdaptor("hsapiens", "GRCh37");
        QueryOptions queryOptions1 = new QueryOptions();

        Query query1 = new Query();
        query1.put(ClinicalDBAdaptor.QueryParams.TRAIT.key(), "alzheimer");
        queryOptions1.add(QueryOptions.INCLUDE, "annotation.traitAssociation.id");
        QueryResult<Variant> queryResult1 = clinicalDBAdaptor.get(query1, queryOptions1);
        assertEquals(1, queryResult1.getNumResults());
        assertTrue(containsAccession(queryResult1, "RCV000172777"));

        Query query2 = new Query();
        query2.put(ClinicalDBAdaptor.QueryParams.TRAIT.key(), "myelofibrosis");
        QueryOptions queryOptions2 = new QueryOptions();
        queryOptions2.add(QueryOptions.INCLUDE, "annotation.traitAssociation.id");
        QueryResult queryResult2 = clinicalDBAdaptor.nativeGet(query2, queryOptions2);
        assertEquals(1, queryResult2.getNumResults());

        Query query4 = new Query();
        query4.put(ClinicalDBAdaptor.QueryParams.REGION.key(),
                new Region("2", 170360030, 170362030));
        QueryOptions queryOptions4 = new QueryOptions();
        queryOptions4.add(QueryOptions.INCLUDE, "annotation.traitAssociation.id");
        QueryResult<Variant> queryResult4 = clinicalDBAdaptor.get(query4, queryOptions4);
        assertEquals(2, queryResult4.getNumTotalResults());
        assertTrue(containsAccession(queryResult4, "COSM4624460"));
        assertTrue(containsAccession(queryResult4, "RCV000171500"));

        Query query5 = new Query();
        query5.put(ClinicalDBAdaptor.QueryParams.CLINICALSIGNIFICANCE.key(), "likely_pathogenic");
        QueryOptions queryOptions5 = new QueryOptions();
        QueryResult queryResult5 = clinicalDBAdaptor.nativeGet(query5, queryOptions5);
        assertEquals(2, queryResult5.getNumTotalResults());

        Query query6 = new Query();
        query6.put(ClinicalDBAdaptor.QueryParams.FEATURE.key(), "APOE");
        QueryOptions queryOptions6 = new QueryOptions();
        queryOptions6.put(QueryOptions.SORT, "chromosome,start");
        queryOptions6.put(QueryOptions.INCLUDE, "chromosome,start,annotation.consequenceTypes.geneName,annotation.traitAssociation.genomicFeatures.xrefs.symbol,annotation.consequenceTypes,annotation.traitAssociation.id");
        QueryResult queryResult6 = clinicalDBAdaptor.nativeGet(query6, queryOptions6);
        // Check sorted output
        int previousStart = -1;
        for (Document document : (List<Document>) queryResult6.getResult()) {
            assertTrue(previousStart < document.getInteger("start"));
        }

        queryOptions6.remove(QueryOptions.SORT);
        query6.put(ClinicalDBAdaptor.QueryParams.SOURCE.key(), "clinvar");
        QueryResult queryResult7 = clinicalDBAdaptor.nativeGet(query6, queryOptions6);
        assertEquals(1, queryResult7.getNumTotalResults());

        query6.put(ClinicalDBAdaptor.QueryParams.SOURCE.key(), "cosmic");
        QueryResult<Variant> queryResult8 = clinicalDBAdaptor.get(query6, queryOptions6);
        assertEquals(1, queryResult8.getNumTotalResults());
        List<String> geneSymbols = queryResult8.getResult().get(0).getAnnotation().getTraitAssociation().stream()
                .map((evidenceEntry) -> evidenceEntry.getGenomicFeatures().get(0).getXrefs().get("symbol"))
                .collect(Collectors.toList());
        geneSymbols.addAll(queryResult8.getResult().get(0).getAnnotation().getConsequenceTypes().stream()
                .map((consequenceType) -> consequenceType.getGeneName())
                .collect(Collectors.toList()));
        assertThat(geneSymbols, CoreMatchers.hasItem("APOE"));

        Query query7 = new Query();
        query7.put(ClinicalDBAdaptor.QueryParams.ACCESSION.key(),"COSM306824");
        query7.put(ClinicalDBAdaptor.QueryParams.SOURCE.key(), "cosmic");
        QueryOptions options = new QueryOptions();
        QueryResult<Variant> queryResult9 = clinicalDBAdaptor.get(query7, options);
        assertNotNull("Should return the queryResult of id=COSM306824", queryResult9.getResult());
        assertThat(queryResult9.getResult().get(0).getAnnotation().getTraitAssociation().stream()
                        .map((evidenceEntry) -> evidenceEntry.getGenomicFeatures().get(0).getXrefs().get("symbol"))
                        .collect(Collectors.toList()),
                CoreMatchers.hasItem("FMN2"));

        // Look for a variant that matches the exact name of a disease
        Query query10 = new Query();
        query10.put(ClinicalDBAdaptor.QueryParams.SOURCE.key(),"clinvar");
        query10.put(ClinicalDBAdaptor.QueryParams.TRAIT_EXACT_MATCH.key(),"hypercholesterolemia and hypertriglyceridemia, type iii");
        QueryOptions options10 = new QueryOptions();
        QueryResult queryResult10 = clinicalDBAdaptor.nativeGet(query10, options10);
        assertEquals(1, queryResult10.getNumResults());


    }

    private boolean containsAccession(QueryResult<Variant> queryResult1, String accession) {
        boolean found = false;
        int i = 0;
        while (i < queryResult1.getNumResults() && !found) {
            int j = 0;
            while (j < queryResult1.getResult().get(i).getAnnotation().getTraitAssociation().size() && !found) {
                found = queryResult1
                        .getResult()
                        .get(i)
                        .getAnnotation()
                        .getTraitAssociation()
                        .get(j)
                        .getId()
                        .equals(accession);
                j++;
            }
            i++;
        }
        return found;
    }

    @Test
    public void getDiseases() throws Exception {
        ClinicalDBAdaptor clinicalDBAdaptor = dbAdaptorFactory.getClinicalDBAdaptor("hsapiens", "GRCh37");

        // Get all diseases
        Query query = new Query();
        QueryOptions queryOptions = new QueryOptions();
        QueryResult queryResult = clinicalDBAdaptor.getDiseases(query,queryOptions);
        assertEquals(7, queryResult.getNumResults());

        // Get diseases checks trait
        Query query2 = new Query();
        query2.put(ClinicalDBAdaptor.QueryParams.SOURCE.key(), "clinvar");
        query2.put(ClinicalDBAdaptor.QueryParams.TRAIT.key(), "disease");
        QueryOptions queryOptions2 = new QueryOptions();
        QueryResult queryResult2 = clinicalDBAdaptor.getDiseases(query2,queryOptions2);
        assertEquals(2, queryResult2.getNumResults());

        Iterator itr = queryResult2.getResult().iterator();
        System.out.println("Test: getDiseases() checks trait 'disease':");
        while (itr.hasNext()) {
            Object element = itr.next();
            System.out.println("\t" + element);
        }

        // Get diseases checks trait with pattern
        Query query3 = new Query();
        query3.put(ClinicalDBAdaptor.QueryParams.SOURCE.key(), "clinvar");
        query3.put(ClinicalDBAdaptor.QueryParams.TRAIT.key(), "hypertriglyceridemia[a-z, ]*type iii");
        QueryOptions queryOptions3 = new QueryOptions();
        QueryResult queryResult3 = clinicalDBAdaptor.getDiseases(query3,queryOptions3);
        assertEquals(1, queryResult3.getNumResults());

        itr = queryResult3.getResult().iterator();
        System.out.println("Test: getDiseases() checks trait and heritable_trait:");
        while (itr.hasNext()) {
            Object element = itr.next();
            System.out.println("\t" + element);
        }
    }
}