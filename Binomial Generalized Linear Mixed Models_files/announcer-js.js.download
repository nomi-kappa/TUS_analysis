/**
  * Announcer WordPress plugin - Main JS file
  * Author : Aakash Chakravarthy ( www.aakashweb.com - aakash.19493@gmail.com )
  *
  * v3.0
  *
  **/

jQuery(document).ready(function(){
	
	var ancr = jQuery( '.announcer' );

	if(ancr.length > 0){
		
		var body = jQuery( 'body' );
		var tfloat = jQuery( '.announcer-top-float' );
		var bfloat = jQuery( '.announcer-bottom-float' );
		var tstatic = jQuery( '.announcer-top-static' );
		
		
		window.ancr_body_pad = {
			'top': parseInt( body.css( 'padding-top' ) ),
			'bottom': parseInt( body.css( 'padding-bottom' ) )
		};
		
		// Wrap it up !!
		if( tstatic.length > 0 ){
			tstatic.prependTo( 'body' );
		}
		
		if( tfloat.length > 0 ){
		
			body.append( '<div id="announcer-top"></div>' );
			tfloat.removeClass( 'announcer-top-float' ).appendTo( '#announcer-top' );

		}
		
		if( bfloat.length > 0 ){
		
			body.append( '<div id="announcer-bottom"></div>' );
			bfloat.removeClass( 'announcer-bottom-float' ).appendTo( '#announcer-bottom' );

		}
		
		ancr.each(function(){
		
			var obj = jQuery( this );
			var id = obj.attr( 'data-id' );
			
			if( announcer_get_cookie( id ) == 'hidden' ){
				
				obj.hide();
				
			}else{
			
				var effect = obj.attr('data-effect');
				var effdur = parseInt(obj.attr('data-effdur'));
				var height = obj.outerHeight();
				var pos = obj.attr('data-pos');
				
				announcer_set_cookie( id, 'visible' );
				
				if( effect == 'slide' ){

					switch( pos ){
						case 'manual': obj.show(); announcer_adjpad(); break;
						default: obj.hide().slideDown( effdur, announcer_adjpad ); break;
					}
					
				}else if( effect == 'fade' ){
					
					obj.hide().fadeIn( effdur, announcer_adjpad );
					
				}
				
			}
		
		});
		
		
		// Clsoe button event
		jQuery('.announcer-closebt').click(function(){
			var id = jQuery(this).parent().attr( 'data-id' );
			
			announcer_set_cookie( id, 'hidden' );

			jQuery(this).parent().slideUp(function(){
				announcer_adjpad();
			});
			
		});
		
		var resizeTimer;
		jQuery(window).resize(function() {
			clearTimeout(resizeTimer);
			resizeTimer = setTimeout(announcer_adjpad, 200);
		});
		
	}
});


function announcer_get_cookie( id ) {
  var value = "; " + document.cookie;
  var parts = value.split("; announcer-" + id + "=");
  if (parts.length == 2) return parts.pop().split(";").shift();
}


function announcer_set_cookie( id, val ){
	var expDate = new Date();
	expDate.setDate(expDate.getDate()+1);
	document.cookie = 'announcer-' + id + '=' + val + '; expires=' + expDate.toGMTString() + '; path=/'
}


function announcer_adjpad(){
	
	var t = jQuery( '#announcer-top' ),
		b = jQuery( '#announcer-bottom' ),
		body = jQuery( 'body' );
		bpad = window.ancr_body_pad;
		
		if( t.length > 0 ) {
			body.animate( { paddingTop: bpad.top + parseInt( t.outerHeight() )  });
		}
		
		if( b.length > 0 ) {
			body.animate( { paddingBottom: bpad.bottom + parseInt( b.outerHeight() ) });
		}
	
}