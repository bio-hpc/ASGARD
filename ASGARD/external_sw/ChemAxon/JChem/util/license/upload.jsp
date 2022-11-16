<%@ page import="java.io.InputStream, chemaxon.util.*,
                 java.io.File,
                 java.io.StringWriter,
                 java.io.PrintWriter,
                 java.util.regex.Pattern,
                 java.util.regex.Matcher,
                 java.io.IOException"%>
 <%--

--%>
<%@ page contentType="text/html;charset=UTF-8" language="java" %>



<html>
  <head><title>Uploading file</title></head>
  <body bgcolor="#e1e1e1">

<%!
static public String stackToString(Throwable e) {
    try {
	StringWriter sw = new StringWriter();
	PrintWriter pw = new PrintWriter(sw);
	e.printStackTrace(pw);
	return sw.toString();
    } catch(Exception e2) {
	return "stackToString error";
    }
}
%>


<%
    String filePath = request.getParameter("filePath");

    Exception e;

    String fileSep=System.getProperty("file.separator");

    String chemDir=DotfileUtil.getDotDir().getAbsolutePath();
    if (chemDir.endsWith(fileSep)) {
        chemDir=chemDir.substring(0, chemDir.length()-1);
    }

    //Creating chemaxon directory if not exists
    File d = new File(chemDir);
    d.mkdirs();

//  Set the file, get parameters from the os.
    String fileName="";
//    String filePath = chemDir +fileSep+ fileName;
    
    Pattern fileNamePattern1=Pattern.compile("(.+)(license.cxl$)");
    Matcher matcher1=fileNamePattern1.matcher(filePath);
    if (matcher1.matches()){
        fileName="license.cxl";
    }
    
    Pattern fileNamePattern2=Pattern.compile("(.+)(licenses.dat$)");
    Matcher matcher2=fileNamePattern2.matcher(filePath);
    if (matcher2.matches()){
        fileName="licenses.dat";
    }
    
    if (fileName == ""){
	e=new IOException();
    }
    
    //New thread to reads the inputstream of the request and saves the file.
    UploadThread uh = new UploadThread();

    uh.setParams(request.getContentType(), request.getContentLength(),
 	request.getInputStream(), chemDir, fileName);

    //Starts the thread, wait until finished.
    uh.start();

    do {
        Thread.sleep(200);
    } while(!uh.isFinished());

    boolean invalid = uh.getFileSize()<=2;

    if(!uh.getSuccess() || invalid) {
%>
    An error occured during the file upload.
<%

        if (invalid) {
        %>
        <p>
        The file does not exist, or its size is zero.

        <p>
        <%
        }

        e=uh.getError();
        if (e!=null) {
%>
            <pre>
Stack trace:
-----------


<%= stackToString(e) %>
            </pre>
<%
        }
    } else {
%>
    File succesfully updated.
    <p>
    <b>NOTE:</b> You have to restart the web server for the license key
        changes to take effect.
<%
    }
%>
    <br>
    <br>
    <input type="button" value="Close this window" onClick="window.close()"/>

  </p>      
  </body>
</html>
