package org.opencb.cellbase.build.transform.serializers.mongodb;

import com.fasterxml.jackson.core.JsonProcessingException;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.databind.ObjectWriter;
import com.fasterxml.jackson.databind.PropertyNamingStrategy;
import com.fasterxml.jackson.databind.SerializationFeature;
import com.fasterxml.jackson.databind.cfg.MapperConfig;
import com.fasterxml.jackson.databind.introspect.AnnotatedField;
import org.opencb.cellbase.build.transform.MutationParser;
import org.opencb.cellbase.build.transform.serializers.CellbaseSerializer;
import org.opencb.cellbase.core.common.GenericFeatureChunk;
import org.opencb.cellbase.core.common.core.Gene;
import org.opencb.cellbase.core.common.core.GenomeSequenceChunk;
import org.opencb.cellbase.core.common.variation.Mutation;
import org.opencb.cellbase.core.common.variation.Variation;

import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: imedina
 * Date: 8/28/13
 * Time: 5:41 PM
 * To change this template use File | Settings | File Templates.
 */
public class MongoDBSerializer implements CellbaseSerializer {

    private File file;
    private Path outdirPath;
    private Path outfilePath;

    private Map<String, BufferedWriter> bufferedWriterMap;

    private BufferedWriter genomeSequenceBufferedWriter;
    // Variation data is too big to be stored in a single file,
    // data is split in different files
    private Map<String, BufferedWriter> variationBufferedWriter;
    private BufferedWriter mutationBufferedWriter;

    private ObjectMapper jsonObjectMapper;
    private ObjectWriter jsonObjectWriter;

    private int chunkSize = 2000;

    public MongoDBSerializer(File file) throws IOException {
        this.file = file;
        init();
    }

    private void init() throws IOException {
        if(file.exists() && file.isDirectory() && file.canWrite()) {
            outdirPath = file.toPath();
        }else {
            outfilePath = file.toPath();
        }

        bufferedWriterMap = new Hashtable<>(50);
        variationBufferedWriter = new HashMap<>(40);

        jsonObjectMapper = new ObjectMapper();
//        jsonObjectMapper.setPropertyNamingStrategy(new GeneNamingStrategy());

        jsonObjectWriter = jsonObjectMapper.writer();

//        PropertyNamingStrategy propertyNamingStrategy = new PropertyNamingStrategy() {
//            @Override
//            public String nameForField(MapperConfig<?> mapperConfig, AnnotatedField annotatedField, String s) {
//                return super.nameForField(mapperConfig, annotatedField, s);    //To change body of overridden methods use File | Settings | File Templates.
//            }
//        };


    }


    @Override
    public void serialize(GenomeSequenceChunk genomeSequenceChunk) {
        try {
            if(genomeSequenceBufferedWriter == null) {
                genomeSequenceBufferedWriter = Files.newBufferedWriter(outfilePath, Charset.defaultCharset());
            }
            genomeSequenceBufferedWriter.write(jsonObjectWriter.writeValueAsString(genomeSequenceChunk));
            genomeSequenceBufferedWriter.newLine();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        };
    }

    @Override
    public void serialize(Gene gene) {
        try {
            if(bufferedWriterMap.get("gene") == null) {
                bufferedWriterMap.put("gene", Files.newBufferedWriter(outdirPath.resolve("gene.json"), Charset.defaultCharset()));
            }
            bufferedWriterMap.get("gene").write(jsonObjectWriter.writeValueAsString(gene));
            bufferedWriterMap.get("gene").newLine();
        } catch (JsonProcessingException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    @Override
    public void serialize(Variation variation) {
        try {
            if(variationBufferedWriter.get(variation.getChromosome()) == null) {
                variationBufferedWriter.put(variation.getChromosome(), Files.newBufferedWriter(outdirPath.resolve("variation_chr" + variation.getChromosome() + ".json"), Charset.defaultCharset()));
            }
            variationBufferedWriter.get(variation.getChromosome()).write(jsonObjectWriter.writeValueAsString(variation));
            variationBufferedWriter.get(variation.getChromosome()).newLine();
        } catch (JsonProcessingException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }

    @Override
    public void serialize(Mutation mutation) {
        try {
            if(mutationBufferedWriter == null) {
                mutationBufferedWriter = Files.newBufferedWriter(outdirPath.resolve("mutation.json"), Charset.defaultCharset());
            }
            mutationBufferedWriter.write(jsonObjectWriter.writeValueAsString(mutation));
            mutationBufferedWriter.newLine();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        };
    }

    @Override
    public void serialize(GenericFeatureChunk genericFeatureChunk) {
        try {
            if(bufferedWriterMap.get("regulatory") == null) {
                bufferedWriterMap.put("regulatory", Files.newBufferedWriter(outdirPath.resolve("regulatory_region.json"), Charset.defaultCharset()));
            }
            bufferedWriterMap.get("regulatory").write(jsonObjectWriter.writeValueAsString(genericFeatureChunk));
            bufferedWriterMap.get("regulatory").newLine();
        } catch (JsonProcessingException e) {
            e.printStackTrace();
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
    }


    @Override
    public void close() {
        String id;
        try {

            closeBufferedWriter(genomeSequenceBufferedWriter);
            closeBufferedWriter(mutationBufferedWriter);
        } catch (IOException e) {
            e.printStackTrace();  //To change body of catch statement use File | Settings | File Templates.
        }
        ;

        Iterator<String> iter = bufferedWriterMap.keySet().iterator();
        while(iter.hasNext()) {
            id = iter.next();
            if(bufferedWriterMap.get(id) != null) {
                try {
                    bufferedWriterMap.get(id).close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

        iter = variationBufferedWriter.keySet().iterator();
        while(iter.hasNext()) {
            id = iter.next();
            if(variationBufferedWriter.get(id) != null) {
                try {
                    variationBufferedWriter.get(id).close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }
    }

    private void closeBufferedWriter(BufferedWriter bufferedWriter) throws IOException {
        if(bufferedWriter != null) {
            bufferedWriter.close();
        }
    }
}