#! /usr/bin/env python

HTML_HEAD = '''<!DOCTYPE html>
<html lang="en">
<head>
  <title> </title>
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <link href="//netdna.bootstrapcdn.com/twitter-bootstrap/2.3.1/css/bootstrap-combined.min.css" rel="stylesheet">
  <style>
    body {
      padding-top: 60px; /* 60px to make the container go all the way to the bottom of the topbar */
    }
  </style>
  <link href="//netdna.bootstrapcdn.com/twitter-bootstrap/2.3.1/css/bootstrap-responsive.min.css" rel="stylesheet">
  <link href="http://synbio.jgi-psf.org:8180/analysis/static/css/pacbio.status.css" rel="stylesheet">
</head>
'''

HTML_BODY = '''<body>
<div class="container-fluid">
  <!-- Page header -->
  <div class="row-fluid">
    <div class="span12">
      <div class="page-header">
        <h1> Analysis: %(aname)s </h1>
      </div>
    </div>
  </div> <!-- end row -->
  
  <div class="row-fluid">
    <!-- Analysis Info -->
    <div class="span6">
      <table class="table table-condensed pretty" id="dtlinks"></table>
    </div> <!-- end span8 -->

<!--
    <div class="span2" style="text-align:center;">
      <a href="./results/%(aname)s.xlsx" class="btn btn-mini" style="margin-top:10px;margin-bottom:20px;">Excel file</a>
      <table class="table table-nolines">
        <tr><th>Perfect</th><td id="npass"></td></tr>
        <tr><th>Fixes</th><td id="nfix"></td></tr>
        <tr><th>Failed</th><td id="nfail"></td></tr>
      </table>
    </div>
-->

    <!-- Legend -->
    <div class="span4 offset8" style="position:fixed;">
      <table class="table table-nolines">
        <thead>
        <tr><th colspan=3>Call </th><tr>
        </thead>
        <tbody>
          <tr><td><span class="badge badge-perfect">F</span></td><td><strong>Flawless</strong></td><td> everything is looking good</td></tr>
          <tr><td><span class="badge badge-almost">A</span></td><td><strong>Almost</strong></td><td> all variants within 10 bases</td></tr>       
          <tr><td><span class="badge badge-incomplete">I</span></td><td><strong>Incomplete</strong></td><td> not all positions covered </td></tr>
          <tr><td><span class="badge badge-lowcov">L</span></td><td><strong>Low coverage</strong></td><td> mean coverage < 30</td></tr>
          <tr><td><span class="badge badge-errors">E</span></td><td><strong>Errors</strong></td><td> variants </td></tr>
          <tr><td><span class="badge badge-dips">D</span></td><td><strong>Dips</strong></td><td> Dips </td></tr>
          <tr><td><span class="badge badge-nocall">?</span></td><td><strong>No call</strong></td><td></td></tr>
        </tbody>
      </table>
    </div> <!-- end span4 -->
    
  </div> <!-- end row -->

  <!-- Results table -->
  <div class="row-fluid">
    <div class="span8">
      <h2> Results </h2>
      <table class="datatable table table-condensed pretty" id="dtresults"></table>
    </div> <!-- end span8 -->
  </div> <!-- end row -->

</div> <!-- end container-fluid -->

</body>
</html>
'''

JS_INCLUDE = '''
  <script src="http://code.jquery.com/jquery.min.js"></script>
  <script src="//netdna.bootstrapcdn.com/twitter-bootstrap/2.3.1/js/bootstrap.min.js"></script>    
  <script type="text/javascript" charset="utf8" src="http://ajax.aspnetcdn.com/ajax/jquery.dataTables/1.9.4/jquery.dataTables.min.js"></script>
  <!-- script type="text/javascript" charset="utf8" src="/static/js/bootstrap-bigmodal.js"></script -->

  <script>
    var analysis_name = "%(aname)s";
    var baseurl = "%(baseurl)s";
    var analysis_url = baseurl + analysis_name + '/';
    var pools = %(pools_json)s;
    var analysis_outcomes = %(outcomes_json)s;
    var references = %(references_json)s;
  </script>
'''

JS_LINK = '''
<script>
  var protocols = %(protocols_json)s;
  var pools = %(pools_json)s;

  var dtLinkCols = [{"sTitle":"Pool","sWidth":"15%%"},{"sTitle":"SampleID","sWidth":"25%%"}];
  $.each(protocols,function(j,protocol) {
    dtLinkCols.push({"sTitle":protocol});
  });
  dtLinkCols.push({"sTitle":"IGV links","sWidth":"25%%"});
  
  var dtLinks = [];
  $.each(pools,function(j,p){
    var dt_row = ['<b>'+p.name+'</b>',p.sample_id];
    $.each(protocols,function(j,protocol) {
      if(protocol in p.jobs) {
        var td_str = '<a href="' + analysis_url + p.name + '/' + p.jobs[protocol] + '/call_summary.txt">' + p.jobs[protocol] + '</a> ';
        dt_row.push(td_str);
      }
      else {
        var td_str = "-";
        dt_row.push(td_str);
      }
    });
    var td_str = '<small><a href="' + analysis_url + 'results/' + p.name + '.igv.xml">igv.xml</a></small>'
    dt_row.push(td_str);
    dtLinks.push(dt_row);
  });
  
  $(function() {
    $("#dtlinks").dataTable({
        'bPaginate':false,
        'bFilter':false,
        'aaSorting':[],
        'bInfo':false,        
        'aoColumns': dtLinkCols, //array of columns with "sTitle"
        'aaData':dtLinks // array of arrays with row data
    });
  });
</script>
'''

"""
  var analysis_outcomes = %(outcomes_json)s;
  var references = %(references_json)s;
  var bestbets = %(bestbets_json)s;
"""

JS_RESULT = '''
<script>
  /*
  var calculateTotals = function() {
    var counts = {'fail':0,'fix':0,'pass':0};
    $.each(bestbets,function(k,v){
      if(!v) {
        ++counts.fail;
      } else {
        if(analysis_outcomes[k][bestbets[k]].call=="perfect") ++counts.pass;
        else if (analysis_outcomes[k][bestbets[k]].call=="almost") ++counts.fix;
      }
    });
    var nrefs = counts.pass + counts.fix + counts.fail;
    $("#npass").html(counts.pass + ' (' + (100*counts.pass/nrefs).toFixed(1) +'%%)');
    $("#nfix").html(counts.fix + ' (' + (100*counts.fix/nrefs).toFixed(1) +'%%)');
    $("#nfail").html(counts.fail + ' (' + (100*counts.fail/nrefs).toFixed(1) +'%%)');
  }
  */
  
  var dtCols = [{"sTitle":"Reference","sWidth":"30%%"},{"sTitle":"Length","sWidth":"10%%"}];
  $.each(pools,function(i,p){
    var hstr = p.name;
    dtCols.push({"sTitle":hstr}); 
  });
  
  var dtData = [];
  $.each(references,function(i,ref) {
    var newrefname = ref.name;
    if (newrefname.length>30){ newrefname = newrefname.substring(0,30) + "...";}
    var reftip = '<span data-toggle="tooltip" title="' + ref.name + '">' + newrefname + '</span>';
    var rowdata = [reftip,ref.length];
    //var bbet = 'pool1'; //bestbets[ref.name];
    $.each(pools,function(j,p){
      //var td_str = '<div class="dropdown"><a class="dropdown-toggle" role="button" data-toggle="dropdown" data-target="#" href=".">'
      var outcome = analysis_outcomes[ref.name][p.name];      
      var cstr = '';
      if(outcome) {
        cstr = '<a href="'+ outcome.igv +'"><span class="badge badge-' + outcome.call + ' badge-digit3">' + outcome.display + '</span></a>';
        //if(p.name==bbet) {
        //  cstr = '<span class="badge badge-' + outcome.call + ' badge-digit3" data-bestpool="true">' + outcome.display + '</span>';
        //}
        //else {
        //  
        //}
      }
      else {
        cstr = '<span class="badge">N</span>'
      }
      rowdata.push(cstr);
    });
    dtData.push(rowdata);
  });

  $(function() {
    $("#dtresults").dataTable({
        'bPaginate':false,
        'bFilter':true,
        'aaSorting':[],
        'aoColumns': dtCols, //array of columns with "sTitle"
        'bAutoWidth': false,
        'aaData':dtData // array of arrays with row data
    });
    //$('span[data-bestpool]').each(function(index){$(this).closest("td").addClass("bestclone")});
    //$.each(bestbets,function(k,v){if(!v) {$('span[title="'+ k +'"]').closest("tr").addClass("failrow");}});
    //calculateTotals();
  });
</script>
'''
