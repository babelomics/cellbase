/*
 * Copyright 2015 OpenCB
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

package org.opencb.cellbase.client.rest;

import org.opencb.biodata.formats.protein.uniprot.v201504jaxb.Entry;
import org.opencb.biodata.models.core.Gene;
import org.opencb.biodata.models.core.Transcript;
import org.opencb.biodata.models.core.TranscriptTfbs;
import org.opencb.biodata.models.variant.Variant;
import org.opencb.cellbase.client.config.ClientConfiguration;
import org.opencb.commons.datastore.core.Query;
import org.opencb.commons.datastore.core.QueryResponse;

import java.io.IOException;
import java.util.Map;
import java.util.List;

/**
 * Created by imedina on 12/05/16.
 */
public class GeneClient extends FeatureClient<Gene> {


    public GeneClient(ClientConfiguration clientConfiguration) {
        super(clientConfiguration);

        this.clazz = Gene.class;

        this.category = "feature";
        this.subcategory = "gene";
    }

    public QueryResponse<String> getBiotypes(Query query) throws IOException {
        return execute("biotype", query, String.class);
    }

    @Deprecated
    public QueryResponse<Gene> list(Query query) throws IOException {
        return execute("list", query, Gene.class);
    }

    public QueryResponse<Transcript> getTranscript(String id, Map<String, Object> params) throws IOException {
        return execute(id, "transcript", params, Transcript.class);
    }

    public QueryResponse<TranscriptTfbs> getTfbs(String id, Map<String, Object> params) throws IOException {
        return execute(id, "tfbs", params, TranscriptTfbs.class);
    }

    public QueryResponse<Variant> getSnp(List<String> id, Map<String, Object> params) throws IOException {
        return execute(id, "snp", params, Variant.class);
    }

    public QueryResponse<Entry> getProtein(String id, Map<String, Object> params) throws IOException {
        return execute(id, "protein", params, Entry.class);
    }

//    public QueryResponse<> getClinical(String id, Map<String, Object> params) throws IOException {
//        return execute(id, "clinical", params, );
//    }

}
