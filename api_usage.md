Requesting useful info via the API:
Basic info:
   >>> r.json()["data"][0]["name"]
   u'0.2-um-passable microorganisms in deep-sea hydrothermal fluid'
   >>> r.json()["data"][0]["id"]
   u'mgm4448226.3'
   >>> r.json()["data"][0]["sequence_type"]
   u'WGS'
   >>> r.json()["data"][0]["mixs"]["seq_method"] 
   u'454'
seq_method is what can contain "assembled"

Raw stats:
   >>> pprint([str(key) for key in r.json()["data"][0]["statistics"]["sequence_stats"].keys() if "raw" in key])
   ['average_gc_ratio_raw',
    'ambig_sequence_count_raw',
    'ambig_char_count_raw',
    'average_length_raw',
    'bp_count_raw',
    'standard_deviation_gc_content_raw',
    'drisee_score_raw',
    'standard_deviation_length_raw',
    'average_ambig_chars_raw',
    'sequence_count_raw',
    'standard_deviation_gc_ratio_raw',
    'average_gc_content_raw']

Searching:
  - Get just the "assembled" metagenomes:
      metagenome?limit=1&order=name&verbosity=stats&seq_method=assembled
