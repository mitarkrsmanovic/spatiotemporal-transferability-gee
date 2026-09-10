//ADDING AND FILTERING S2 COLLECTION

var s2Filtered = s2
.filter(ee.Filter.date("2025-05-01", "2025-09-01"))
.filter(ee.Filter.bounds(roi))
.filter(ee.Filter.lte("CLOUDY_PIXEL_PERCENTAGE", 10));


var csPlus = ee.ImageCollection('GOOGLE/CLOUD_SCORE_PLUS/V1/S2_HARMONIZED');

var QA_BAND = 'cs_cdf';

var CLEAR_THRESHOLD = 0.60;

var composite = s2
    .filterBounds(roi)
    .filterDate('2025-05-01', '2025-10-01')
    .linkCollection(csPlus, [QA_BAND])
    .map(function(img) {
      return img.updateMask(img.select(QA_BAND).gte(CLEAR_THRESHOLD));
    })
    .median()
    .multiply(0.0001);

var rgbVis = { 
min: 0,
max: 3000,
bands: ["B4", "B3", "B2"]
};


//CREATIN POLYGONS (PALE, SOKOLAC, ROGATICA)

var palePolygon = pale_polygon.geometry()


var polygonVis = { 
color: "FF0000",
fillColor: "00000000",
width: 2
};

var paleFiltered = gaul
.filter(ee.Filter.eq("ADM1_NAME", "Republika Srpska"))
.filter(ee.Filter.eq("ADM2_NAME", "Pale"));

var paleFilteredGeometry = paleFiltered.geometry();

var paleArea = paleFilteredGeometry.area();
var paleAreaKm2 = paleArea.divide(1e6);
//print(paleAreaKm2);


//print(s2Filtered.size());
//Map.addLayer(paleFiltered.style(polygonVis2), {}, "Pale Polygon");
//Map.addLayer(composite, rgbVis, "RGB composite");
//Map.addLayer(pale_polygon.style(polygonVis), { }, "Pale Polygon");
//Map.centerObject(palePolygon, 11);


var sokolacFiltered = gaul
.filter(ee.Filter.eq("ADM1_NAME", "Republika Srpska"))
.filter(ee.Filter.eq("ADM2_NAME", "Sokolac"));

var sokolacPolygon = sokolacFiltered.geometry();

var sokolacArea = sokolacPolygon.area();
var sokolacAreaKm2 = sokolacArea.divide(1e6);
//print(sokolacAreaKm2);

var polygonVis2 = { 
color: "FF0000",
fillColor: "00000000",
width: 2
};

//Map.addLayer(sokolacFiltered.style(polygonVis2), {}, "Sokolac Polygon");


var rogaticaFiltered = gaul
.filter(ee.Filter.eq("ADM1_NAME", "Republika Srpska"))
.filter(ee.Filter.eq("ADM2_NAME", "Rogatica"));

var rogaticaPolygon = rogaticaFiltered.geometry();

var polygonVis3 = { 
color: "FF0000",
fillColor: "00000000",
width: 2
};

//Map.addLayer(rogaticaFiltered.style(polygonVis3), {}, "Rogatica Polygon");



//CALCULATE INDICES (NDBI, MNDWI, NDRE1, SAVI, EVI)

var indices = function(image) {

var mndwi = image.normalizedDifference(["B3", "B11"]).rename("mndwi");
var ndre1 = image.normalizedDifference(["B8A", "B5"]).rename("ndre1");
var urbanindex = image.expression(
"(SWIR2 - NIR) / (SWIR2 + NIR)", {
"SWIR2": image.select("B12"),
"NIR": image.select("B8")
}).rename("urban_index");
var bsi = image.expression(
"((SWIR2 + RED) - (NIR + BLUE)) / ((SWIR2 + RED) + (NIR + BLUE))", { 
"SWIR2": image.select("B12"),
"RED": image.select("B4"),
"NIR": image.select("B8"),
"BLUE": image.select("B2")}).rename("bsi");
var savi = image.expression(
"((NIR - RED) / (NIR + RED +L )) * (1 + L)", { 
"NIR": image.select("B8"),
"RED": image.select("B4"),
"L": 0.5}).rename("savi");
var evi = image.expression(
"2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))", { 
"NIR": image.select("B8"),
"RED": image.select("B4"),
"BLUE": image.select("B2")}).rename("evi");

return image.addBands([mndwi, ndre1, urbanindex, bsi, savi, evi]);
};


//ADDING INDICES TO COMPOSITE IMAGE

var compositeIndices = indices(composite).select(["B2","B3", "B4", "B8","B11", "B12", "mndwi", "ndre1", "urban_index","bsi", "savi", "evi"]);


//ADDING AND FILTERING S1 COLLECTION

var s1Filtered = s1
.filter(ee.Filter.date("2025-05-01", "2025-09-01"))
.filter(ee.Filter.bounds(roi))
.filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VV"))
.filter(ee.Filter.listContains("transmitterReceiverPolarisation", "VH"))
.filter(ee.Filter.eq("instrumentMode", "IW"));

var desc = s1Filtered
.filter(ee.Filter.eq("orbitProperties_pass", "DESCENDING"));

var asc = s1Filtered
.filter(ee.Filter.eq("orbitProperties_pass", "ASCENDING"));

var descVV = desc.select("VV").mean().rename("VV_desc");
var descVH = desc.select("VH").mean().rename("VH_desc");
var ascVV = asc.select("VV").mean().rename("VV_asc");
var ascVH = asc.select("VH").mean().rename("VH_asc");

var s1Composite = ee.Image.cat([descVV, descVH, ascVV, ascVH]);

var vvMean = s1Composite.expression("(VV_desc + VV_asc) / 2", {
  "VV_desc": s1Composite.select("VV_desc"),
  "VV_asc": s1Composite.select("VV_asc")
}).rename("VV");

var vhMean = s1Composite.expression("(VH_desc + VH_asc) / 2", {
  "VH_desc": s1Composite.select("VH_desc"),
  "VH_asc": s1Composite.select("VH_asc")
}).rename("VH");

var ratio = vvMean.divide(vhMean).rename("VV_VH_ratio");

var sarFinal = ee.Image.cat([vvMean, vhMean, ratio]);


//ADDING DEM (ELEVATION + SLOPE)

var slope = ee.Terrain.slope(dem).rename("slope");


// ADDING TEXTURE

var b8_smooth = composite.select("B8").focalMean(10, 'circle', 'meters');
var gray = b8_smooth.multiply(255).toInt16();
var glcm = gray.glcmTexture({size: 7});
var textureBands = glcm.select(["B8_asm", "B8_contrast", "B8_savg"])
                       .rename(["asm", "contrast", "savg"]);


//Adding Sesonal Composites

var springComposite = s2
.filter(ee.Filter.date("2025-03-01", "2025-06-01"))
.filter(ee.Filter.bounds(roi))
.filter(ee.Filter.lte("CLOUDY_PIXEL_PERCENTAGE", 10))

var springCompositeFC = springComposite
.linkCollection(csPlus, [QA_BAND])
.map(function(img) { 
return img.updateMask(img.select(QA_BAND).gte(CLEAR_THRESHOLD));
})
.median()
.multiply(0.0001);
//Map.addLayer(springCompositeFC, rgbVis, "Spring Composite");

var springNDVI = springCompositeFC.normalizedDifference(["B8", "B4"]).rename("spring_ndvi");
var springBSI = springCompositeFC.expression(
"((SWIR2 + RED) - (NIR + BLUE)) / ((SWIR2 + RED) + (NIR + BLUE))", { 
"SWIR2": springCompositeFC.select("B12"),
"RED": springCompositeFC.select("B4"),
"NIR": springCompositeFC.select("B8"),
"BLUE": springCompositeFC.select("B2")}).rename("spring_bsi");
var springEVI = springCompositeFC.expression(
"2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))", { 
"NIR": springCompositeFC.select("B8"),
"RED": springCompositeFC.select("B4"),
"BLUE": springCompositeFC.select("B2")}).rename("spring_evi");

var autumComposite = s2
.filter(ee.Filter.date("2025-09-01", "2025-11-01"))
.filter(ee.Filter.bounds(roi))
.filter(ee.Filter.lte("CLOUDY_PIXEL_PERCENTAGE", 10))

var autumCompositeFC = autumComposite
.linkCollection(csPlus, [QA_BAND])
.map(function(img) { 
return img.updateMask(img.select(QA_BAND).gte(CLEAR_THRESHOLD));
})
.median()
.multiply(0.0001);
//Map.addLayer(autumCompositeFC, rgbVis, "Autum Composite");

var autumNDVI = autumCompositeFC.normalizedDifference(["B8", "B4"]).rename("autum_ndvi");
var autumBSI = autumCompositeFC.expression(
"((SWIR2 + RED) - (NIR + BLUE)) / ((SWIR2 + RED) + (NIR + BLUE))", { 
"SWIR2": autumCompositeFC.select("B12"),
"RED": autumCompositeFC.select("B4"),
"NIR": autumCompositeFC.select("B8"),
"BLUE": autumCompositeFC.select("B2")}).rename("autum_bsi");
var autumEVI = autumCompositeFC.expression(
"2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))", { 
"NIR": autumCompositeFC.select("B8"),
"RED": autumCompositeFC.select("B4"),
"BLUE": autumCompositeFC.select("B2")}).rename("autum_evi");


//Creating Final Samples Image (Composite + Indices + DEM + Spatial Context)

var finalBands = compositeIndices.addBands(dem.rename("elevation")).addBands(slope).addBands(springNDVI).addBands(springBSI).addBands(autumNDVI).addBands(autumBSI).addBands(textureBands).addBands(springEVI).addBands(autumEVI).addBands(sarFinal);


//CREATING SAMPLE COLLECTION

var sampleCollection = forest.merge(grassland).merge(build_up_land).merge(bare_rock).merge(crops).merge(bare_land).merge(sparse_vegetation);


//CREATING SAMPLING GRID

var researchArea = ee.FeatureCollection([
ee.Feature(palePolygon),
ee.Feature(sokolacPolygon),
ee.Feature(rogaticaPolygon)
]).union();

//Map.addLayer(researchArea, {color: "red"}, "Research Area");


var removeAutocorrelation = function(inputSamples, gridSize) { 

var proj = ee.Projection("EPSG:32634");
var gridImage = ee.Image.pixelCoordinates(proj).reproject(proj, null, gridSize);
var pointsWithGrid = gridImage.reduceRegions({
collection: inputSamples, 
reducer: ee.Reducer.first(), 
scale: gridSize
});

var uniquePoints = pointsWithGrid.map(function(f) {
var x = ee.Number(f.get('x')).format('%.0f');
var y = ee.Number(f.get('y')).format('%.0f');
var gridId = ee.String(x).cat('_').cat(ee.String(y));
return f.set('grid_id', gridId);
});

return uniquePoints.distinct('grid_id');
};

var cleanSampleCollection = removeAutocorrelation(sampleCollection, 50);

//print("PALE - Originalni broj uzoraka:", sampleCollection.size());
//print("PALE - Očišćen broj uzoraka (50m razmak):", cleanSampleCollection.size());




//Creating Training Samples and Validation Samples Collections

//var randomSamples = sampleCollection.randomColumn();
//var trainingData = randomSamples.filter(ee.Filter.lt("random", 0.6));
//var validationData = randomSamples.filter(ee.Filter.gte("random", 0.6));


var classes = ee.List([1,2,3,4,5,6,7]);

// Kreiramo trening i validaciju ODVOJENO za svaku klasu unutar mape
var trainingData = ee.FeatureCollection(classes.map(function(c) {
return cleanSampleCollection
.filter(ee.Filter.eq("landcover", c))
.randomColumn("random", 42)
.filter(ee.Filter.lt("random", 0.6));
})).flatten();

var validationData = ee.FeatureCollection(classes.map(function(c) {
return cleanSampleCollection
.filter(ee.Filter.eq("landcover", c))
.randomColumn("random", 42)
.filter(ee.Filter.gte("random", 0.6));
})).flatten();


//print("Trening po klasama:", trainingData.aggregate_histogram("landcover"));
//print("Validacija po klasama:", validationData.aggregate_histogram("landcover"));


//Classifying LULC Image with RF Model

var trainingSamples = finalBands.sampleRegions({
collection: trainingData, 
properties: ["landcover"], 
scale: 10,
geometries: true
});

var predictionBands = ["B2","B3", "B4", "B11", "savi", "evi", "urban_index", "slope", "autum_bsi", "VV", "VH", "VV_VH_ratio"];

var classifier = ee.Classifier.smileRandomForest({numberOfTrees: 90,  minLeafPopulation: 1, bagFraction: 0.8}).train({
features: trainingSamples, 
classProperty: "landcover", 
inputProperties: predictionBands});

var classifiedImage = finalBands.classify(classifier);

var classesVis = { 
min: 1,
max: 7,
palette: ["#1a9850", "#91cf60", "#d73027", "#f5f5f5", "#fe9929", "#cc4c02", "#99d8c9"]
};

//Map.addLayer(classifiedImage.clip(polygon), classesVis, "lulc 2025");


//Export Pale LULC Map

//Export.image.toDrive({
//image: classifiedImage.clip(polygon).uint8(), 
//description: "pale_25_lulc", 
//folder: "earthengine", 
//fileNamePrefix: "pale_lulc_25", 
//region: polygon, 
//scale: 10, 
//maxPixels: 1e13});


//Creating Validation Samples

var validationSamples = finalBands.sampleRegions({
collection: validationData, 
properties: ["landcover"], 
scale: 10});


//Classifying Validation Samples

var validationDataClassified = validationSamples.classify(classifier, 'class_rf_pale');


//Error Matrix Validation(Confusion Matrix, Accuracy, Kappa, Consumers Accuracy, Producers Accuracy)

var testConfusionMatrix = validationDataClassified.errorMatrix("landcover", "class_rf_pale");
print("Pale RF Confusion Matrix", testConfusionMatrix);
//print("Pale RF Accuracy", testConfusionMatrix.accuracy());
//print("Pale RF Kappa", testConfusionMatrix.kappa());
//print("Pale RF Consumers Accuracy", testConfusionMatrix.consumersAccuracy());
//print("Pale RF Producers Accuracy", testConfusionMatrix.producersAccuracy());


//Exporting Error Matrix Results

var confusionMatrixExport = ee.Dictionary({ 
"Overall Accuracy": testConfusionMatrix.accuracy(),
"Kappa": testConfusionMatrix.kappa(),
"Consumers Accuracy": testConfusionMatrix.consumersAccuracy(),
"Producers Accuracy": testConfusionMatrix.producersAccuracy(),
'Matrix': testConfusionMatrix.array()
});

var tcmFeature = ee.FeatureCollection([
ee.Feature(null, confusionMatrixExport)
]);

//Export.table.toDrive({
//collection: tcmFeature, 
//description: "Error_Matrix", 
//folder: "earthengine", 
//fileNamePrefix: "Pale_LULC_Accuracy_Report_2025", 
//fileFormat: "CSV"
//});


//HYPER PARAMETERS TUNING

var numberOfTreesList = ee.List.sequence(10, 150, 10);
var bagFractionList = ee.List.sequence(0.1, 0.9, 0.1);
var minLeafList = ee.List([1, 5, 10]);

var accuracies = numberOfTreesList.map(function(numTrees) {
return bagFractionList.map(function(bagFraction) {
return minLeafList.map(function(minLeaf) {
var classifier = ee.Classifier.smileRandomForest({
numberOfTrees: numTrees,
bagFraction: bagFraction,
minLeafPopulation: minLeaf
}).train({
features: trainingSamples, 
classProperty: "landcover", 
inputProperties: predictionBands
});

var accuracy = validationSamples
.classify(classifier)
.errorMatrix("landcover", "classification")
.accuracy();

return ee.Feature(null, {
        "accuracy": accuracy,
        "numberOfTrees": numTrees,
        "bagFraction": bagFraction,
        "minLeafPopulation": minLeaf
      });
    });
  });
}).flatten().flatten();

var result = ee.FeatureCollection(accuracies);
var sortedResult = result.sort("accuracy", false);
var bestModel = result.sort("accuracy", false).first();
//print("Best Model", bestModel);


//Exporting Hyperparameters Tuning Results

//Eksportovanje tabele na Google Drive
//Export.table.toDrive({
//collection: sortedResult,
//description: 'RF_Hyperparameter_Tuning_Results',
//fileFormat: 'CSV',
//selectors: ['accuracy', 'numberOfTrees', 'bagFraction', 'minLeafPopulation'] // Određuje redoslijed kolona u CSV-u
//});


//Finding Misclassified points

//var misClassified = validationDataClassified.filter(
//ee.Filter.neq("landcover", "classification")
//);
//print("Total Misclassified Points", misClassified.size());
//print("Details of Misclassified Points", misClassified);

//Map.addLayer(misClassified, {color: "FF00FF"}, "Misclassified Points");

//var errorList = misClassified.map(function(feature) {
//return feature.select(['landcover', 'classification']);
//});
//print("Sample of Error Labels (Truth vs Predicted):", errorList.limit(10));

//Map.addLayer(misClassified, {color: "blue"}, "Misclassified Points (154)");

//var errorAnalysisMatrix = misClassified.errorMatrix("landcover", "classification");
//print("Summary of Misclassifications (Truth vs. Prediction):", errorAnalysisMatrix);


//Creating Cluster Image

var seedsPaleRF = ee.Algorithms.Image.Segmentation.seedGrid(4);

var snicPaleRF = ee.Algorithms.Image.Segmentation.SNIC({
image: finalBands.select("B.*"), 
size: 4, 
compactness: 0, 
connectivity: 4, 
neighborhoodSize: 10, 
seeds: seedsPaleRF
});

var clustersPaleRF = snicPaleRF.select("clusters");

var smoothedPaleRF = classifiedImage.addBands(clustersPaleRF);

var clustersMajority = smoothedPaleRF.reduceConnectedComponents({
reducer: ee.Reducer.mode(),
labelBand: "clusters"
});

//Map.addLayer(clustersMajority.clip(polygon), classesVis, "Processed Using Clusters");


//Export.image.toDrive({
//image: clustersMajority.clip(palePolygon).uint8(),
//description: "pale_RF_25_lulc", 
//folder: "earthengine", 
//fileNamePrefix: "pale_RF_lulc_25", 
//region: palePolygon, 
//scale: 10,
//maxPixels: 1e13});


//Calculating Area(km2) of ROI

var polygonArea = palePolygon.area();
//print("Polygon Area", polygonArea.divide(1e6));


//Calculating Area(km2) of LULC Classes

var areaImagePaleRF = ee.Image.pixelArea().divide(1e6).addBands(clustersMajority);
//print("Area Pale RF Image", areaImagePaleRF);

var areaPaleRF = areaImagePaleRF.reduceRegion({
reducer: ee.Reducer.sum().group({
groupField: 1, 
groupName: "classification"}),
geometry: palePolygon,
scale: 10,
maxPixels: 1e10, 
tileScale: 8
});
//print(areaPaleRF);


//Exporting Areas Data

var list = ee.List(areaPaleRF.get("groups"));
//print("List", list);

var featureCollection = ee.FeatureCollection(list.map(function(item) { 
var dict = ee.Dictionary(item);

return ee.Feature(null, { 
"Class_id": dict.get("classification"),
"Area_km2": dict.get("sum")
});
}));

//Export.table.toDrive({
//collection: featureCollection, 
//description: "Area_km2_2025", 
//folder: "earthengine", 
//fileNamePrefix: "Area_Pale_RF_2025"});



//FINAL VALIDATION OF CLASSIFICATION RESULTS

var validationDataSokolac = forest_2.merge(grassland_2).merge(build_up_land_2).merge(bare_rock_2).merge(crops_2).merge(bare_land_2).merge(sparse_vegetation_2);

var cleanSokolacSamples = removeAutocorrelation(validationDataSokolac, 50);

//print("SOKOLAC - Originalni broj uzoraka:", validationDataSokolac.size());
//print("SOKOLAC - Očišćen broj uzoraka (50m razmak):", cleanSokolacSamples.size());

var sampledDataSokolac = finalBands.sampleRegions({
  collection: cleanSokolacSamples,
  properties: ["landcover"],
  scale: 10
});

var sokolacRFModel = sampledDataSokolac.classify(classifier, 'class_rf_sokolac')     

var classifiedImageSokolac = finalBands.classify(classifier);

var seedsSokolacRF = ee.Algorithms.Image.Segmentation.seedGrid(4);

var snicSokolacRF = ee.Algorithms.Image.Segmentation.SNIC({
image: finalBands.select("B.*"), 
size: 4, 
compactness: 0, 
connectivity: 4, 
neighborhoodSize: 10, 
seeds: seedsSokolacRF
});

var clustersSokolacRF = snicSokolacRF.select("clusters");

var smoothedSokolac = classifiedImageSokolac.addBands(clustersSokolacRF);

var clustersMajoritySokolac = smoothedSokolac.reduceConnectedComponents({
reducer: ee.Reducer.mode(), 
labelBand: "clusters"
});


//Export.image.toDrive({
//image: clustersMajoritySokolac.clip(sokolacPolygon).uint8(), 
//description: "sokolac_RF_25_lulc", 
//folder: "earthengine", 
//fileNamePrefix: "sokolac_RF_lulc_25", 
//region: sokolacPolygon, 
//scale: 10,
//maxPixels: 1e13});


var testConfusionMatrixSokolac = sokolacRFModel.errorMatrix("landcover", "class_rf_sokolac");
print("Sokolac RF Confusion Matrix", testConfusionMatrixSokolac);
//print("Sokolac RF Accuracy", testConfusionMatrixSokolac.accuracy());
//print("Sokolac RF Kappa", testConfusionMatrixSokolac.kappa());
//print("Sokolac RF Consumers Accuracy", testConfusionMatrixSokolac.consumersAccuracy());
//print("Sokolac RF Producers Accuracy", testConfusionMatrixSokolac.producersAccuracy());

var testConfusionMatrixSokolacExport = ee.Dictionary({ 
"Overall Accuracy": testConfusionMatrixSokolac.accuracy(),
"Kappa": testConfusionMatrixSokolac.kappa(),
"Consumers Accuracy": testConfusionMatrixSokolac.consumersAccuracy(),
"Producers Accuracy": testConfusionMatrixSokolac.producersAccuracy(),
'Matrix': testConfusionMatrixSokolac.array()
});

var tcmFeature = ee.FeatureCollection([
ee.Feature(null, testConfusionMatrixSokolacExport)
]);

//Export.table.toDrive({
//collection: tcmFeature, 
//description: "Error_Matrix", 
//folder: "earthengine", 
//fileNamePrefix: "Sokolac_LULC_Accuracy_Report_2025", 
//fileFormat: "CSV"
//});


//Calculating Area(km2) of LULC Classes

var areaImageSokolacRF = ee.Image.pixelArea().divide(1e6).addBands(clustersMajoritySokolac);
//print("Area Image Sokolac", areaImageSokolacRF);

var areaSokolacRF = areaImageSokolacRF.reduceRegion({
reducer: ee.Reducer.sum().group({
groupField: 1, 
groupName: "classification"}),
geometry: sokolacPolygon,
scale: 10,
maxPixels: 1e10, 
tileScale: 8
});
//print(areaSokolacRF);


//Exporting Areas Data

var list = ee.List(areaSokolacRF.get("groups"));
//print("List", list);

var featureCollection = ee.FeatureCollection(list.map(function(item) { 
var dict = ee.Dictionary(item);

return ee.Feature(null, { 
"Class_id": dict.get("classification"),
"Area_km2": dict.get("sum")
});
}));

//Export.table.toDrive({
//collection: featureCollection, 
//description: "Area_km2_2025", 
//folder: "earthengine", 
//fileNamePrefix: "Area_Sokolac_RF_2025"});


// ============================================================
// RF FEATURE IMPORTANCE
// ============================================================

var rfExplain = classifier.explain();
var rfImportance = ee.Dictionary(rfExplain.get('importance'));

var rfSum = rfImportance.values().reduce(ee.Reducer.sum());
var rfImportancePct = rfImportance.map(function(key, val) {
  return ee.Number(val).divide(rfSum).multiply(100);
});

var rfBandNames = rfImportancePct.keys();

var rfImportanceList = rfBandNames.map(function(key) {
  return ee.Feature(null, {
    'band': key,
    'importance_pct': rfImportancePct.get(key)
  });
});

var rfImportanceTable = ee.FeatureCollection(rfImportanceList);

// Eksportuj u Google Drive
//Export.table.toDrive({
//collection: rfImportanceTable.sort('importance_pct', false),
//description: 'RF_Feature_Importance_Final',
//fileFormat: 'CSV',
//selectors: ['band', 'importance_pct']
//});



//GBM CLASSIFIER

var gbmClassifier = ee.Classifier.smileGradientTreeBoost({
numberOfTrees: 60, 
shrinkage: 0.1, 
samplingRate: 0.8, 
maxNodes: null, 
loss: "LeastAbsoluteDeviation"
}).train({
features: trainingSamples, 
classProperty:"landcover" , 
inputProperties: predictionBands});

var classifiedImage2 = finalBands.classify(gbmClassifier);

var classesVis = { 
min: 1,
max: 7,
palette: ["#1a9850", "#91cf60", "#d73027", "#f5f5f5", "#fe9929", "#cc4c02", "#99d8c9"]
};

//Map.addLayer(classifiedImage2.clip(polygon), classesVis, "lulc 2025");

//Creating Validation Samples

var gbmValidationSamples = finalBands.sampleRegions({
collection: validationData, 
properties: ["landcover"], 
scale: 10});


//Classifying Validation Samples

var gbmValidationDataClassified = gbmValidationSamples.classify(gbmClassifier, 'class_gbm_pale');


//Creating Cluster Image

var seedsPaleGBM = ee.Algorithms.Image.Segmentation.seedGrid(4);

var snicPaleGBM = ee.Algorithms.Image.Segmentation.SNIC({
image: finalBands.select("B.*"), 
size: 4, 
compactness: 0, 
connectivity: 4, 
neighborhoodSize: 10, 
seeds: seedsPaleGBM
});

var clustersPaleGBM = snicPaleGBM.select("clusters");

var smoothedPaleGBM = classifiedImage2.addBands(clustersPaleGBM);

var gbmClustersMajority = smoothedPaleGBM.reduceConnectedComponents({
reducer: ee.Reducer.mode(), 
labelBand: "clusters"
});

//Map.addLayer(gbmClustersMajority.clip(polygon), classesVis, "Processed Using Clusters");


//Export.image.toDrive({
//image: gbmClustersMajority.clip(palePolygon).uint8(), 
//description: "pale_GBM_25_lulc", 
//folder: "earthengine", 
//fileNamePrefix: "pale_GBM_lulc_25",
//region: palePolygon, 
//scale: 10,
//maxPixels: 1e13});


//Calculating Area(km2) of LULC Classes

var areaImagePaleGBM = ee.Image.pixelArea().divide(1e6).addBands(gbmClustersMajority);
//print("Area Pale RF Image", areaImagePaleRF);

var areaPaleGBM = areaImagePaleGBM.reduceRegion({
reducer: ee.Reducer.sum().group({
groupField: 1, 
groupName: "classification"}),
geometry: palePolygon,
scale: 10,
maxPixels: 1e10, 
tileScale: 8
});
//print(areaPaleRF);


//Exporting Areas Data

var list = ee.List(areaPaleGBM.get("groups"));
//print("List", list);

var featureCollection = ee.FeatureCollection(list.map(function(item) { 
var dict = ee.Dictionary(item);

return ee.Feature(null, { 
"Class_id": dict.get("classification"),
"Area_km2": dict.get("sum")
});
}));

//Export.table.toDrive({
//collection: featureCollection, 
//description: "Area_km2_2025", 
//folder: "earthengine", 
//fileNamePrefix: "Area_Pale_GBM_2025"});


//Error Matrix Validation(Confusion Matrix, Accuracy, Kappa, Consumers Accuracy, Producers Accuracy)

var gbmTestConfusionMatrix = gbmValidationDataClassified.errorMatrix("landcover", "class_gbm_pale");
print("Pale GBM Confusion Matrix", gbmTestConfusionMatrix);
//print("Pale GBM Accuracy", gbmTestConfusionMatrix.accuracy());
//print("Pale GBM Kappa", gbmTestConfusionMatrix.kappa());
//print("Pale GBM Consumers Accuracy", gbmTestConfusionMatrix.consumersAccuracy());
//print("Pale GBM Producers Accuracy", gbmTestConfusionMatrix.producersAccuracy());


//Hyper Parameters Tuning

var numberOfTreesList = ee.List.sequence(10, 150, 10);
var samplingRateList = ee.List([0.5, 0.8]);
var shrinkageList = ee.List([0.01, 0.1]);

var accuracies = numberOfTreesList.map(function(numTrees) {
return samplingRateList.map(function(sampleRate) {
return shrinkageList.map(function(shrink) {
var classifier = ee.Classifier.smileGradientTreeBoost({
numberOfTrees: numTrees,
samplingRate: sampleRate,
shrinkage: shrink
}).train
({
features: trainingSamples, 
classProperty: "landcover", 
inputProperties: predictionBands
});

var accuracy = gbmValidationSamples
.classify(classifier)
.errorMatrix("landcover", "classification")
.accuracy();

return ee.Feature(null, {
"accuracy": accuracy,
"numberOfTrees": numTrees,
"samplingRate": sampleRate,
"shrinkage": shrink
});
});
});
}).flatten().flatten();

var result = ee.FeatureCollection(accuracies);
var sortedResult = result.sort("accuracy", false);
var bestModel = result.sort("accuracy", false).first();
//print("Best Model", bestModel);


//Eksportovanje u CSV format na Google Drive
//Export.table.toDrive({
//collection: sortedResult,
//description: 'GBM_Hyperparameter_Tuning_Results',
//fileFormat: 'CSV',
//selectors: ['accuracy', 'numberOfTrees', 'samplingRate', 'shrinkage']
//});


//Creating Cluster Image

var gbmClassifiedImageSokolac = finalBands.classify(gbmClassifier);

var gbmSeedsSokolac = ee.Algorithms.Image.Segmentation.seedGrid(4);

var gbmSnicSokolac = ee.Algorithms.Image.Segmentation.SNIC({
image: finalBands.select("B.*"), 
size: 4, 
compactness: 0, 
connectivity: 4, 
neighborhoodSize: 10, 
seeds: gbmSeedsSokolac
});

var gbmClustersSokolac = gbmSnicSokolac.select("clusters");

var gbmSmoothedSokolac = gbmClassifiedImageSokolac.addBands(gbmClustersSokolac);

var gbmClustersMajoritySokolac = gbmSmoothedSokolac.reduceConnectedComponents({
reducer: ee.Reducer.mode(), 
labelBand: "clusters"
});

//Map.addLayer(gbmClustersMajority.clip(polygon), classesVis, "Processed Using Clusters");


//Export.image.toDrive({
//image: gbmClustersMajoritySokolac.clip(sokolacPolygon).uint8(), 
//description: "sokolac_GBM_25_lulc", 
//folder: "earthengine", 
//fileNamePrefix: "sokolac_GBM_lulc_25", 
//region: sokolacPolygon, 
//scale: 10,
//maxPixels: 1e13});

var sokolacGBMModel = sampledDataSokolac.classify(gbmClassifier, 'class_gbm_sokolac')     

var classifiedImageSokolac = finalBands.classify(gbmClassifier);


//Calculating Area(km2) of LULC Classes

var areaImageSokolacGBM = ee.Image.pixelArea().divide(1e6).addBands(gbmClustersMajoritySokolac);
//print("Area Image Sokolac", areaImageSokolacRF);

var areaSokolacGBM = areaImageSokolacGBM.reduceRegion({
reducer: ee.Reducer.sum().group({
groupField: 1, 
groupName: "classification"}),
geometry: sokolacPolygon,
scale: 10,
maxPixels: 1e10, 
tileScale: 8
});
//print(areaSokolacGBM);


//Exporting Areas Data

var list = ee.List(areaSokolacGBM.get("groups"));
//print("List", list);

var featureCollection = ee.FeatureCollection(list.map(function(item) { 
var dict = ee.Dictionary(item);

return ee.Feature(null, { 
"Class_id": dict.get("classification"),
"Area_km2": dict.get("sum")
});
}));

//Export.table.toDrive({
//collection: featureCollection, 
//description: "Area_km2_2025", 
//folder: "earthengine", 
//fileNamePrefix: "Area_Sokolac_GBM_2025"});

var gbmTestConfusionMatrixSokolac = sokolacGBMModel.errorMatrix("landcover", "class_gbm_sokolac");
print("Sokolac GBM Confusion Matrix", gbmTestConfusionMatrixSokolac);
//print("Sokolac GBM Accuracy", gbmTestConfusionMatrixSokolac.accuracy());
//print("Sokolac GBM Kappa", gbmTestConfusionMatrixSokolac.kappa());
//print("Sokolac GBM Consumers Accuracy", gbmTestConfusionMatrixSokolac.consumersAccuracy());
//print("Sokolac GBM Producers Accuracy", gbmTestConfusionMatrixSokolac.producersAccuracy());



// 1. Funkcija za pakovanje metrika u jedan Feature
var createAccuracyFeature = function(confusionMatrix, areaName) {
return ee.Feature(null, {
'Area': areaName,
'Overall_Accuracy': confusionMatrix.accuracy(),
'Kappa': confusionMatrix.kappa(),
'Producers_Accuracy': confusionMatrix.producersAccuracy().toList().flatten(),
'Consumers_Accuracy': confusionMatrix.consumersAccuracy().toList().flatten(),
'Matrix': confusionMatrix.array()
});
};

// 2. Kreiranje feature-a za Pale i Sokolac
var paleMetrics = createAccuracyFeature(gbmTestConfusionMatrix, "Pale");
var sokolacMetrics = createAccuracyFeature(gbmTestConfusionMatrixSokolac, "Sokolac");

// 3. Spajanje u FeatureCollection
var accuracyReport = ee.FeatureCollection([paleMetrics, sokolacMetrics]);

// 4. Eksportovanje u Google Drive (CSV format)
//Export.table.toDrive({
//collection: accuracyReport,
//description: 'GBM_LULC_Accuracy_Report_2025',
//folder: 'earthengine',
//fileNamePrefix: 'LULC_GBM_Metrics_Pale_Sokolac_2025',
//fileFormat: 'CSV',
//selectors: ['Area', 'Overall_Accuracy', 'Kappa', 'Producers_Accuracy', 'Consumers_Accuracy']
//});


var gbmExplain = gbmClassifier.explain();
var gbmImportance = ee.Dictionary(gbmExplain.get('importance'));

var gbmSum = gbmImportance.values().reduce(ee.Reducer.sum());
var importancePct = gbmImportance.map(function(key, val) {
  return ee.Number(val).divide(gbmSum).multiply(100);
});

var bandNames = importancePct.keys();

var importanceList = bandNames.map(function(key) {
  return ee.Feature(null, {
    'band': key,
    'importance_pct': importancePct.get(key)
  });
});

var importanceTable = ee.FeatureCollection(importanceList);


Export.table.toDrive({
collection: importanceTable.sort('importance_pct', false),
description: 'GBM_Feature_Importance_Final',
fileFormat: 'CSV',
selectors: ['band', 'importance_pct']
});


// ============================================================
// SEPARABILITY ANALYSIS – sve klase (1–7)
// ============================================================

var separabilityBands = ["B2","B3","B4","B11","savi","evi",
                         "urban_index","slope","autum_bsi",
                         "VV","VH","VV_VH_ratio"];

var MIN_SAMPLES = 30;

// Nazivi klasa za interpretaciju u outputu
var classNames = {
  1: "forest",
  2: "grassland",
  3: "build_up_land",
  4: "bare_rock",
  5: "crops",
  6: "bare_land",
  7: "sparse_vegetation"
};

// ============================================================
// 1. UZORKOVANJE PODATAKA
// ============================================================
var paleSamples = finalBands.select(separabilityBands)
  .sampleRegions({
    collection: trainingData,
    properties: ["landcover"],
    scale: 10,
    geometries: true
  });

var sokolacSamples = finalBands.select(separabilityBands)
  .sampleRegions({
    collection: cleanSokolacSamples,
    properties: ["landcover"],
    scale: 10,
    geometries: true
  });

// ============================================================
// 2. DOMAIN SHIFT ANALIZA – sve klase
// ============================================================
var dsClasses = [1, 2, 3, 4, 5, 6, 7];

var domainShiftResults = dsClasses.map(function(c) {
  var paleCls    = paleSamples.filter(ee.Filter.eq("landcover", c));
  var sokolacCls = sokolacSamples.filter(ee.Filter.eq("landcover", c));

  var paleCount    = paleCls.size();
  var sokolacCount = sokolacCls.size();

  var paleStats = separabilityBands.map(function(band) {
    return ee.Feature(null, {
      "klasa":      c,
      "class_name": classNames[c],
      "band":       band,
      "mean":       paleCls.aggregate_mean(band),
      "stdDev":     paleCls.aggregate_total_sd(band),
      "domain":     "pale",
      "n":          paleCount
    });
  });

  var sokolacStats = separabilityBands.map(function(band) {
    return ee.Feature(null, {
      "klasa":      c,
      "class_name": classNames[c],
      "band":       band,
      "mean":       sokolacCls.aggregate_mean(band),
      "stdDev":     sokolacCls.aggregate_total_sd(band),
      "domain":     "sokolac",
      "n":          sokolacCount
    });
  });

  return ee.FeatureCollection(paleStats.concat(sokolacStats));
});

var domainShiftFC = ee.FeatureCollection(domainShiftResults).flatten();
//print("Domain Shift – statistika po bandu i domenu:", domainShiftFC);

// ============================================================
// 3. EFFECT SIZE – Cohen's d za sve klase
// ============================================================
var effectSizeResults = dsClasses.map(function(c) {
  var paleCls    = paleSamples.filter(ee.Filter.eq("landcover", c));
  var sokolacCls = sokolacSamples.filter(ee.Filter.eq("landcover", c));

  var n1 = paleCls.size();
  var n2 = sokolacCls.size();

  var enoughSamples = ee.Number(n1).gte(MIN_SAMPLES)
                        .and(ee.Number(n2).gte(MIN_SAMPLES));

  var cohenDPerBand = separabilityBands.map(function(band) {

    var m1 = ee.Number(paleCls.aggregate_mean(band));
    var m2 = ee.Number(sokolacCls.aggregate_mean(band));
    var s1 = ee.Number(paleCls.aggregate_total_sd(band));
    var s2 = ee.Number(sokolacCls.aggregate_total_sd(band));

    // Standardna pooled SD
    var n1m1      = ee.Number(n1).subtract(1);
    var n2m1      = ee.Number(n2).subtract(1);
    var denom     = ee.Number(n1).add(n2).subtract(2);
    var pooledVar = n1m1.multiply(s1.pow(2))
                        .add(n2m1.multiply(s2.pow(2)))
                        .divide(denom);
    var pooledStd = pooledVar.sqrt();

    var cohenD = ee.Algorithms.If(
      pooledStd.gt(0).and(enoughSamples),
      m1.subtract(m2).abs().divide(pooledStd),
      ee.Number(-9999)
    );

    var dNum = ee.Number(cohenD);
    var interpretation = ee.Algorithms.If(
      dNum.lt(0),         "INVALID – premalo uzoraka",
      ee.Algorithms.If(
        dNum.lt(0.2),     "zanemarljiv",
        ee.Algorithms.If(
          dNum.lt(0.5),   "mali",
          ee.Algorithms.If(
            dNum.lt(0.8), "srednji",
                          "VELIKI – domain shift!"
    ))));

    return ee.Feature(null, {
      "klasa":          c,
      "class_name":     classNames[c],
      "band":           band,
      "cohen_d":        cohenD,
      "pale_mean":      m1,
      "sokolac_mean":   m2,
      "pale_std":       s1,
      "sokolac_std":    s2,
      "pale_n":         n1,
      "sokolac_n":      n2,
      "interpretation": interpretation
    });
  });

  return ee.FeatureCollection(cohenDPerBand);
});

var effectSizeFC = ee.FeatureCollection(effectSizeResults).flatten();
//print("Effect Size (Cohen's d) – po klasi i bandu:", effectSizeFC);

// ============================================================
// 4. PROVJERA UZORAKA – sve klase
// ============================================================
dsClasses.forEach(function(c) {
  var pN = paleSamples.filter(ee.Filter.eq("landcover", c)).size();
  var sN = sokolacSamples.filter(ee.Filter.eq("landcover", c)).size();
  //print("Klasa " + c + " (" + classNames[c] + ") – Pale n =", pN);
  //print("Klasa " + c + " (" + classNames[c] + ") – Sokolac n =", sN);
});

// ============================================================
// 5. EXPORT
// ============================================================
Export.table.toDrive({
collection:  effectSizeFC,
description: "CohenD_DomainShift_sve_klase_2025",
fileFormat:  "CSV",
folder:      "earthengine",
selectors:   ["klasa","class_name","band","cohen_d","interpretation",
"pale_mean","sokolac_mean", "pale_std","sokolac_std", "pale_n","sokolac_n"]});

Export.table.toDrive({
collection: domainShiftFC,
description: "DomainShift_Statistika_sve_klase_2025",
fileFormat: "CSV",
folder: "earthengine",
selectors: ["klasa","class_name","band","domain","mean","stdDev","n"]
});


//EXPORT CLASSIFIERS FOR 2025.

//Export.classifier.toAsset({
//classifier: classifier, 
//description: "RF_Model_Pale_2025", 
//assetId: "users/mitarkrsmanovic/RF_Model_Pale_2025"
//});

//Export.classifier.toAsset({
//classifier: gbmClassifier, 
//description: "GBM_Model_Pale_2025", 
//assetId: "users/mitarkrsmanovic/GBM_Model_Pale_2025"
//});


//var trainingSamplesFiltered = trainingSamples.filter(ee.Filter.notNull(['.geo']));

//Export.table.toAsset({
//collection: trainingSamples,
//description: "TrainingSamples_Pale_2025_Scaled",
//assetId: "users/mitarkrsmanovic/TrainingSamples_Pale_2025_Scaled"
//});



// EXPORT FINAL SAMPLES AND POLYGONS


var paleWithLabel = cleanSampleCollection.map(function(f) {
return f.set('site', 'Pale');
});

var sokolacWithLabel = cleanSokolacSamples.map(function(f) {
return f.set('site', 'Sokolac');
});

// Spoji obe kolekcije
var allSamples = paleWithLabel.merge(sokolacWithLabel);

// Eksportuj kao jedan shapefile
//Export.table.toDrive({
//collection: allSamples,
//description: 'Samples_Pale_Sokolac_2025',
//folder: 'earthengine',
//fileNamePrefix: 'Samples_Pale_Sokolac_2025',
//fileFormat: 'SHP'
//});

// Pale – sve geometrije pretvorene u Polygon
var palePoly = ee.Feature(
paleFiltered.geometry().buffer(1).dissolve(),
{name: 'Pale'}
);

// Sokolac – sve geometrije pretvorene u Polygon
var sokolacPoly = ee.Feature(
sokolacFiltered.geometry().buffer(1).dissolve(),
{name: 'Sokolac'}
);

// Eksportuj zajedno
//Export.table.toDrive({
//collection: ee.FeatureCollection([palePoly, sokolacPoly]),
//description: 'Polygons_Pale_Sokolac_Final',
//folder: 'earthengine',
//fileNamePrefix: 'Polygons_Pale_Sokolac',
//fileFormat: 'SHP'
//});
