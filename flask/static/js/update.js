/**
 * Created by rsw on 4/18/18.
 */
function update_request(){
    updateFilterseq();
    updateStatus();
    $.ajax({
        global: false,
        type: 'POST',
        url: 'update_request', // missing quotes
        dataType: 'json',
        data: JSON.stringify({
            id: $("#reqid").val(),
            priority: $("#priority").val(),
            type: 'priority'
        }),
        contentType: "application/json;charset=utf-8",
        success: function (result) {
            console.log(result);
        },
        error: function (request, status, error) {
            serviceError();
        }
    });
    location.reload()
};

function updateFilterseq(){
    var filters = ["ifu", "r", "i", "g", "u"];
    var arrayLength = filters.length;
    var flt_list = []
    var exptime_list = []

    for (var i = 0; i < arrayLength; i++) {
        if (document.getElementById(filters[i]+'_checked').checked) {
            var repeat = parseInt(document.getElementById(filters[i]+'_repeat').value);
            var exptime = parseInt(document.getElementById(filters[i]+'_exptime').value);
            if ( repeat >= 1 && exptime >= 1)  {
                flt_list.push(repeat+filters[i]);
                exptime_list.push(exptime);
            }

        }

    }
    var flt_seq = '{'+flt_list.join(", ")+'}';
    var exp_seq = '{'+exptime_list.join(", ")+'}';
    $.ajax({
        global: false,
        type: 'POST',
        url: 'update_request', // missing quotes
        dataType: 'json',
        data: JSON.stringify({
            id: $("#reqid").val(),
            filter_seq: flt_seq,
            exptime_seq: exp_seq,
            type: 'filters',
        }),
        contentType: "application/json;charset=utf-8",
        success: function (result) {
            console.log(result);
        },
        error: function (request, status, error) {
            serviceError();
        }
    });
}

function updateStatus() {
     var radioChecked = 'False';
     var checkedStatus = '';
     if (document.getElementById('completed').checked) {
        checkedStatus = 'COMPLETED'
     }
     else if (document.getElementById('observed').checked){
         checkedStatus = 'OBSERVED'
     }
     else if (document.getElementById('pending').checked){
         checkedStatus = 'PENDING'
     }
     else {
         return
     }
     alert(checkedStatus)
     $.ajax({
            global: false,
            type: 'POST',
            url: 'update_request', // missing quotes
            dataType: 'json',
            data: JSON.stringify({
                id: document.getElementById('reqid').value,
                status: checkedStatus,
                type: 'status',
            }),
            contentType: "application/json;charset=utf-8",
            success: function (result) {
                console.log(result);
            },
            error: function (request, status, error) {
                serviceError();
            }
        });
}
