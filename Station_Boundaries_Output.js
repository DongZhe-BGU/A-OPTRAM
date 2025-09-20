var POINTS_RAW = roi;
var WC         = ee.Image('ESA/WorldCover/v200/2021').select('Map');
var SCALE      = 10;
var BUFFER_M   = 1000;
var ID_FIELD   = 'folder';
var MIN_AREA_M2 = 0;
var EXPORT_NAME = '';
var INCLUDE_CIRCLE_BUFFER = false;

function safeName(s) { return String(s).replace(/[^\w\-]+/g, '_'); }

var POINTS = POINTS_RAW
  .map(function(f){ return ee.Feature(f).set('geom_type', f.geometry().type()); })
  .filter(ee.Filter.eq('geom_type', 'Point'))
  .map(function(f){
    var pt  = f.geometry().centroid(1);
    var lon = ee.Number(pt.coordinates().get(0));
    var lat = ee.Number(pt.coordinates().get(1));
    var key = lon.format('%.6f').cat('_').cat(lat.format('%.6f'));
    return ee.Feature(f).set('xy_key', key);
  })
  .distinct('xy_key');

print('POINTS.size (after clean & dedup) =', POINTS.size());

function buildOne(feat) {
  var pt  = feat.geometry().centroid(1);
  var buf = pt.buffer(BUFFER_M);

  var hasIdField = ee.List(feat.propertyNames()).contains(ID_FIELD);
  var lon = ee.Number(pt.coordinates().get(0));
  var lat = ee.Number(pt.coordinates().get(1));
  var sid = ee.String(ee.Algorithms.If(
    hasIdField, feat.getString(ID_FIELD),
    ee.String('pt_').cat(lon.format('%.5f')).cat('_').cat(lat.format('%.5f'))
  ));

  var sampled   = WC.sample(pt, SCALE).first();
  var lcAtPoint = ee.Algorithms.If(sampled, ee.Feature(sampled).get('Map'), null);
  var noLC      = ee.Algorithms.IsEqual(lcAtPoint, null);

  var sameClass = ee.Image(ee.Algorithms.If(
    noLC, ee.Image(0), WC.eq(ee.Number(lcAtPoint)).selfMask()
  )).clip(buf);

  var pixCount = sameClass.reduceRegion({
    reducer: ee.Reducer.count(), geometry: buf, scale: SCALE, bestEffort: true, maxPixels: 1e9
  }).get('Map');
  var hasPix = ee.Number(ee.Algorithms.If(noLC, 0, pixCount)).gt(0);

  var emptyFC = ee.FeatureCollection([]);
  var circleFC = ee.FeatureCollection([
    ee.Feature(buf).copyProperties(feat).set({
      station_id: sid,
      buffer_m  : BUFFER_M,
      area_m2   : buf.area(1),
      area_ha   : buf.area(1).divide(1e4),
      geom_type : 'circle_buffer'
    })
  ]);

  var outFC = ee.FeatureCollection(ee.Algorithms.If(
    hasPix,
    (function() {
      var vectors = sameClass.reduceToVectors({
        geometry: buf,
        scale: SCALE,
        geometryType: 'polygon',
        eightConnected: true,
        bestEffort: true,
        maxPixels: 1e9
      });

      if (MIN_AREA_M2 > 0) {
        vectors = vectors.map(function(f){
          var g = ee.Feature(f).geometry();
          return ee.Feature(g).set('area_raw', g.area(1));
        }).filter(ee.Filter.greaterThanOrEquals('area_raw', MIN_AREA_M2));
      }

      var vSize = vectors.size();
      return ee.Algorithms.If(
        vSize.gt(0),
        (function() {

          var withA = vectors.map(function(f){
            var g = ee.Feature(f).geometry();
            return ee.Feature(g).set('area_raw', g.area(1));
          });
          var biggest = ee.Feature(withA.sort('area_raw', false).first());
          var g = biggest.geometry();
          var a = g.area(1);

          return ee.Algorithms.If(
            a.gt(0),
            ee.FeatureCollection([
              ee.Feature(g)
                .copyProperties(feat)
                .set({
                  station_id : sid,
                  lc_code    : lcAtPoint,
                  buffer_m   : BUFFER_M,
                  area_m2    : a,
                  area_ha    : a.divide(1e4),
                  geom_type  : 'wc_sameclass_maxpatch'
                })
            ]),
            ee.Algorithms.If(INCLUDE_CIRCLE_BUFFER, circleFC, emptyFC)
          );
        })(),
        ee.Algorithms.If(INCLUDE_CIRCLE_BUFFER, circleFC, emptyFC)
      );
    })(),
    ee.Algorithms.If(INCLUDE_CIRCLE_BUFFER, circleFC, emptyFC)
  ));

  return outFC;
}

var perPoint = POINTS.map(buildOne).flatten();
print('Generated features (all) =', perPoint.size());

var perPoint_valid = perPoint
  .map(function(f){
    var g = ee.Feature(f).geometry();
    return ee.Feature(f).set({
      area_ok: g.area(1),
      gtype  : g.type()
    });
  })
  .filter(ee.Filter.gt('area_ok', 0))
  .filter(ee.Filter.inList('gtype', ['Polygon', 'MultiPolygon']));

print('Features after final filter =', perPoint_valid.size());

Map.centerObject(POINTS, 7);
Map.addLayer(perPoint_valid, {}, 'Same-class polygons (valid)');

var SELECTORS = ['station_id', 'lc_code', 'buffer_m', 'area_m2', 'area_ha', 'geom_type'];
var RENAME_TO = ['station',   'lc_code', 'buffer_m', 'area_m2', 'area_ha', 'geom_type'];

var perPoint_out = perPoint_valid.map(function(f){
  var props = ee.Dictionary({
    'station'  : f.get('station_id'),
    'lc_code'  : f.get('lc_code'),
    'buffer_m' : ee.Number(f.get('buffer_m')),
    'area_m2'  : ee.Number(f.get('area_m2')),
    'area_ha'  : ee.Number(f.get('area_ha')),
    'geom_type': f.get('geom_type')
  });
  var g = ee.Feature(f).geometry();
  return ee.Feature(g, props);
});

print('Props check (first feature):', ee.Feature(perPoint_out.first()).propertyNames());

Export.table.toDrive({
  collection    : perPoint_out,
  description   : EXPORT_NAME + '_shp',
  fileNamePrefix: safeName(EXPORT_NAME) + '_shp',
  fileFormat    : 'SHP'
});

print('Exports prepared. Start them in the Tasks panel.');
