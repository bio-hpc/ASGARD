<%--
  JSP file for setting licenses for the Web Server
--%>
<%@ page contentType="text/html;charset=UTF-8" language="java" %>
<html>
<head><title>Setting license for Web Server</title></head>

<body bgcolor="#e1e1e1">

<script LANGUAGE="JavaScript1.1">
<!--
function upload() {
    var form=document.frm;
    var fileName=form.filePath.value;
    if (fileName=="") {
        alert("Specify the new license file.")
        return false;
    }
    else {
        var regExp = new RegExp("(.+)(license.cxl$)|(licenses.dat$)");
        if (regExp.test(filePath.value) == false){
    	    alert("We can accept ChemAxon license files. (license.cxl or licenses.dat)")
            return false;
        }
    }
    
    conf=confirm("Warning: if the license file exists on the server then it will be overwritten.")
    if (!conf) {
        return false;
    }
    form.action="upload.jsp?filePath="+filePath.value;
    form.submit();
    
    return true;
}
-->
</script>

<center><h2>Setting license file for the web server</h2></center>

<br>
<br>
<br>

<form name="frm" ENCTYPE="multipart/form-data" action="upload.jsp" method="post">

<table cellpadding="10" >
<tr>

<td>New license file:</td>
<td><input type="file" size="80" name="filePath" id="filePath"></td>
</tr>
</table>

<input type=button value="Upload file" name="filePath" id="filePath" onClick="upload()">
</form>    
</body>
</html>
