var SHP_FC        = ee.FeatureCollection(roi);
var STATION_FIELD = 'station';
var OUT_PREFIX    = '';

var EXPORT_SELECTORS_BASE = [
  'year','month',
  'NDVI','SAVI_L050','STR'
];

var YEARS  = [];
var MONTHS = [];

var cloudPct = 10;
var scale = 30;
var numPixelsPerMonth = 80000;

var USE_MASK_LC         = true;
var FILTER_NDVI_NONNEG  = true;
var STR_MIN_SWIR        = 0.02;
var STR_CLAMP_MIN       = 0.0;
var STR_CLAMP_MAX       = 30.0;

function maskS2clouds(img) {
  var scl = img.select('SCL');
  var clear = scl.neq(3)
               .and(scl.neq(8))
               .and(scl.neq(9))
               .and(scl.neq(10))
               .and(scl.neq(0))
               .and(scl.neq(1))
               .and(scl.neq(2))
               .and(scl.neq(7))
               .and(scl.neq(11));
  return img.updateMask(clear).copyProperties(img, ['system:time_start']);
}

function toReflectance(img) {
  var bands = ['B2','B3','B4','B5','B6','B7','B8','B8A','B11','B12'];
  var scaled = img.select(bands).toFloat().divide(10000);
  return img.addBands(scaled, null, true).copyProperties(img, img.propertyNames());
}

function saviWithL(nir, red, L) {
  return nir.subtract(red)
            .multiply(1 + L)
            .divide(nir.add(red).add(L));
}

function addIndices(img) {
  var nir  = img.select('B8');
  var red  = img.select('B4');
  var swir = img.select('B12');

  var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');

  var savi_l050 = saviWithL(nir, red, 0.5).rename('SAVI_L050');
  var savi_l025 = saviWithL(nir, red, 0.25).rename('SAVI_L025');
  var savi_l100 = saviWithL(nir, red, 1.0).rename('SAVI_L100');

  var twoNirPlus1 = nir.multiply(2).add(1);
  var radicand = twoNirPlus1.pow(2).subtract(nir.subtract(red).multiply(8));
  radicand = radicand.max(0);
  var msavi = twoNirPlus1.subtract(radicand.sqrt()).divide(2).rename('MSAVI');

  var swirValid = swir.gt(STR_MIN_SWIR).and(swir.lt(1));
  var strRaw = img.expression('(1 - SWIR) * (1 - SWIR) / (2 * SWIR)', {SWIR: swir});
  var str = strRaw.updateMask(swirValid)
                  .clamp(STR_CLAMP_MIN, STR_CLAMP_MAX)
                  .rename('STR');

  return img.addBands([ndvi, savi_l050, savi_l025, savi_l100, msavi, str]);
}

var nonWaterNonUrbanMask = ee.Image(1);
if (USE_MASK_LC) {
  var landCover = ee.Image('ESA/WorldCover/v200/2021').select('Map');
  nonWaterNonUrbanMask = landCover.neq(80).and(landCover.neq(50)).rename('LCmask');
}

var startYear = ee.Number(ee.List(YEARS).reduce(ee.Reducer.min()));
var endYear   = ee.Number(ee.List(YEARS).reduce(ee.Reducer.max()));
var ymList = ee.List(YEARS).map(function(y){
  y = ee.Number(y);
  return ee.List(MONTHS).map(function(m){
    return ee.Dictionary({year: y, month: ee.Number(m)});
  });
}).flatten();

function buildSamplesForAOI(AOI, stationStr) {
  var baseCol = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filterBounds(AOI)
    .filterDate(ee.Date.fromYMD(startYear, 1, 1),
                ee.Date.fromYMD(endYear.add(1), 1, 1))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloudPct))
    .map(maskS2clouds)
    .map(toReflectance)
    .map(addIndices)
    .map(function(img){
      var d = ee.Date(img.get('system:time_start'));
      return img.set({year: d.get('year'), month: d.get('month')});
    })
    .filter(ee.Filter.inList('year', YEARS))
    .filter(ee.Filter.inList('month', MONTHS))
    .map(function(img){ return img.updateMask(nonWaterNonUrbanMask); });

  var allResults = ee.FeatureCollection(ymList.map(function(dm){
    dm = ee.Dictionary(dm);
    var y = ee.Number(dm.get('year'));
    var m = ee.Number(dm.get('month'));

    var start = ee.Date.fromYMD(y, m, 1);
    var end   = start.advance(1, 'month');

    var colMonth = baseCol.filterDate(start, end);
    var count = colMonth.size();
    var monthly = colMonth.median();

    return ee.Algorithms.If(
      count.eq(0),
      ee.FeatureCollection([]),
      (function(){
        var withLL = monthly.addBands(ee.Image.pixelLonLat());

        var toSample = withLL.select([
          'NDVI','SAVI_L050','SAVI_L025','SAVI_L100','MSAVI','STR',
          'longitude','latitude'
        ]);

        var samp = toSample.sample({
          region: AOI,
          scale: scale,
          numPixels: numPixelsPerMonth,
          seed: y.add(m),
          geometries: false
        });

        var feats = samp.map(function(f){
          return ee.Feature(null, {
            station   : stationStr,
            year      : y,
            month     : m,
            longitude : f.get('longitude'),
            latitude  : f.get('latitude'),
            NDVI      : f.get('NDVI'),
            SAVI_L050 : f.get('SAVI_L050'),
            SAVI_L025 : f.get('SAVI_L025'),
            SAVI_L100 : f.get('SAVI_L100'),
            MSAVI     : f.get('MSAVI'),
            STR       : f.get('STR'),
            n_img     : count
          });
        });

        return ee.FeatureCollection(
          FILTER_NDVI_NONNEG ? feats.filter(ee.Filter.gte('NDVI', 0)) : feats
        );
      })()
    );
  })).flatten();

  return allResults;
}

var EXPORT_SELECTORS = ['station'].concat(EXPORT_SELECTORS_BASE);

var stationList = SHP_FC.aggregate_array(STATION_FIELD).distinct();

stationList.evaluate(function(stations) {
  if (!stations || stations.length === 0) {
    print('No stations found for field: ' + STATION_FIELD);
    return;
  }

  stations.forEach(function(st) {
    var fcStation = SHP_FC.filter(ee.Filter.eq(STATION_FIELD, st));

    fcStation.toList(fcStation.size()).evaluate(function(polyList) {
      if (!polyList || polyList.length === 0) return;

      polyList.forEach(function(poly, idx) {
        var geom   = ee.Feature(poly).geometry();
        var label  = String(st).replace(/[^A-Za-z0-9_\-]+/g, '_');
        var desc   = OUT_PREFIX + '_' + label + '_' + idx;

        var samples = buildSamplesForAOI(geom, st);

        Export.table.toDrive({
          collection: samples,
          description: desc,
          fileFormat: 'CSV',
          selectors: EXPORT_SELECTORS
        });
      });

      print('Scheduled exports for station:', st, 'polygons:', polyList.length);
    });
  });
});

var firstGeom = ee.Feature(SHP_FC.first()).geometry();
var previewCol = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(firstGeom)
  .filterDate(ee.Date.fromYMD(startYear, 1, 1), ee.Date.fromYMD(endYear.add(1), 1, 1))
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloudPct))
  .map(maskS2clouds)
  .map(toReflectance)
  .map(addIndices);

var ndviVis = {min: 0, max: 1, palette: ['blue','white','green']};
Map.centerObject(firstGeom, 10);
Map.addLayer(previewCol.median().select('NDVI').clip(firstGeom), ndviVis, 'NDVI median (first polygon)');
